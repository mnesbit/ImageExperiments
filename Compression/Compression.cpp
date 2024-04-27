#include <iostream>
#include <vector>
#include <set>
#include <format>
#include <algorithm>
#include <functional>
#include <filesystem>
#include <random>
#include <cstdlib>
#include <fstream>
#include <map>
#include "../CompressionLib/inc/MatchingPursuit.h"
#include "../CompressionLib/inc/BitBuffer.h"
#include "../CompressionLib/inc/BasisSet.h"
#include "../CompressionLib/inc/Huffman.h"
#include "../CompressionLib/inc/CompressedImage.h"
#include "../ImageHelper/inc/imgloader.h"
#include "../SimpleMatrix/inc/covariance.h"
#include "../ImageHelper/inc/misc.h"

using namespace img;
using namespace matching;
using namespace compressed;
using namespace bitbuffer;
using namespace basis;

// created using bit allocation strategy outlined in:
// Thakur, V. S., Thakur, K., Gupta, S., & Rao, K. R. (2020). Image?independent optimal non?negative integer bit allocation technique for the DCT?based image transform coders. IET Image Processing, 14(1), 11-24.
// Coefficients then derived from max possible range (i.e. 255 * 8 = 2040 for 8 x 8 patches on 8-bit YUV data)
// adjusted for the observed rate of range contraction e.g. 0.9^coefficient index
// divided by 2^bits allocated from above procedure.
// The DC component quantization is capped at 8.0 and other coefficients are capped at 1.0 divisor
static double s_quantY[32] = {
8,
1,
1,
3,
3,
5,
9,
8,
8,
13,
12,
11,
19,
17,
16,
14,
26,
23,
21,
19,
34,
31,
28,
25,
23,
21,
37,
33,
30,
27,
25,
44,
};

static double s_quantU[32] = {
8,
8,
52,
92,
165,
148,
265,
237,
213,
381,
341,
306,
548,
491,
440,
395,
354,
317,
284,
255,
228,
205,
184,
165,
148,
132,
119,
106,
95,
86,
77,
69,
};

static double s_quantV[32] = {
8,
8,
51,
183,
164,
292,
523,
467,
836,
748,
669,
598,
535,
479,
428,
383,
343,
306,
274,
245,
219,
196,
176,
157,
141,
126,
113,
101,
90,
81,
72,
65,
};

std::vector<std::string> readFiles(std::string path) {
	struct stat sb;
	std::vector<std::string> results;

	for (const auto& entry : std::filesystem::directory_iterator(path)) {

		std::filesystem::path outfilename = entry.path();
		std::string outfilename_str = outfilename.string();
		const char* path = outfilename_str.c_str();

		if (stat(path, &sb) == 0) {
			if (!(sb.st_mode & S_IFDIR)) {
				results.push_back(outfilename_str);
			}
			else if ((sb.st_mode & S_IFDIR) == S_IFDIR) {
				const auto subdirectory = readFiles(path);
				results.insert(results.end(), subdirectory.cbegin(), subdirectory.cend());
			}
		}
	}
	return results;
}

void writeDictionaryToPNG(const math::Matrix& dictionary, std::string fileName) {
	size_t blockSize = static_cast<size_t>(sqrt(dictionary.Columns()));
	if (dictionary.Columns() != blockSize * blockSize) {
		throw std::invalid_argument("Dictionary should have square Column count");
	}
	size_t approxWidth = static_cast<size_t>(sqrt(dictionary.Rows()));
	size_t approxHeight = (dictionary.Rows() / approxWidth) + ((dictionary.Rows() % approxWidth == 0) ? 0 : 1);
	image<double>* basisPic = new image<double>(approxWidth * blockSize, approxHeight * blockSize);
	for (size_t i = 0; i < dictionary.Rows(); ++i) {
		size_t blockx = i % approxWidth;
		size_t blocky = i / approxWidth;
		for (size_t x = 0; x < blockSize; ++x) {
			for (size_t y = 0; y < blockSize; ++y) {
				imRef(basisPic, x + blockx * blockSize, y + blocky * blockSize) = dictionary[i][x + (y * blockSize)];
			}
		}
	}
	SaveImageGeneric(basisPic, fileName.c_str(), imgFormat::PNG);
	delete basisPic;
}

int main(int argc, char* argv[]) {
	if (argc != 3) {
		std::cout << "usage: Compression.exe <input file> <output file>" << std::endl;
		return -1;
	}
	const int BlockSize = 8;
	const int K = 32;
	std::cout << "Calculating basis vectors" << std::endl;
	std::vector<Line> basisSet = distinctLineShapes(BlockSize);
	math::Matrix baseDict = createSegmentDictionary(BlockSize, basisSet);
	size_t fullRowCount = baseDict.Rows();
	std::vector<math::Matrix> detailBasisY;
	std::vector<math::Matrix> detailBasisU;
	std::vector<math::Matrix> detailBasisV;
	for (const Line& basis : basisSet) {
		detailBasisY.push_back(createIntraSegmentDictionary(BlockSize, basis, covariancePredictionModelY));
		detailBasisU.push_back(createIntraSegmentDictionary(BlockSize, basis, covariancePredictionModelU));
		detailBasisV.push_back(createIntraSegmentDictionary(BlockSize, basis, covariancePredictionModelV));
	}

	auto dynamic = [&](int prevCoeffs, const std::vector<BasisChoice>& prevChoices, const std::vector<math::Matrix>& detailBasis)  -> math::Matrix {
		size_t atomCount = baseDict.Rows();
		int choice = 0;
		for (int i = 0; i < prevCoeffs; ++i) {
			if (i > 0) {
				choice = choice + zigzagDecode(static_cast<int>(prevChoices[i].deltaId));
			} else {
				choice = static_cast<int>(prevChoices[0].deltaId);
			}
			if (choice >= 0
				&& choice < baseDict.Rows()) {
				atomCount += detailBasis[choice].Rows();
			}
		}
		math::Matrix dictionary(atomCount, BlockSize * BlockSize);
		double* pDict = dictionary.Data();
		memcpy(pDict, baseDict.Data(), sizeof(double) * baseDict.Rows() * baseDict.Columns());
		size_t offset = baseDict.Rows();
		choice = 0;
		for (int i = 0; i < prevCoeffs; ++i) {
			if (i > 0) {
				choice = choice + zigzagDecode(static_cast<int>(prevChoices[i].deltaId));
			} else {
				choice = static_cast<int>(prevChoices[0].deltaId);
			}
			if (choice >= 0
				&& choice < baseDict.Rows()) {
				const math::Matrix& detail = detailBasis[choice];
				memcpy(pDict + (offset * BlockSize * BlockSize), detail.Data(), sizeof(double) * detail.Rows() * detail.Columns());
				offset += static_cast<int>(detail.Rows());
			}
		}
		return dictionary;
		};
	auto dynamicY = [&](int prevCoeffs, const std::vector<BasisChoice>& prevChoices)  -> math::Matrix {
		return dynamic(prevCoeffs, prevChoices, detailBasisY);
		};
	auto dynamicU = [&](int prevCoeffs, const std::vector<BasisChoice>& prevChoices)  -> math::Matrix {
		return dynamic(prevCoeffs, prevChoices, detailBasisU);
		};
	auto dynamicV = [&](int prevCoeffs, const std::vector<BasisChoice>& prevChoices)  -> math::Matrix {
		return dynamic(prevCoeffs, prevChoices, detailBasisV);
		};

	image<rgb>* imgIn = LoadImageGenericRGB(argv[1]);
	std::cout << std::format("start processing image {} using largest {} coefficients.", argv[1], K) << std::endl;
	size_t outputBytesSize;
	std::unique_ptr<uint8_t[]> encodedBytes = encodeImage(
		imgIn,
		K,
		BlockSize,
		s_quantY, s_quantU, s_quantV,
		dynamicY, dynamicU, dynamicV,
		outputBytesSize);

	double squaredError = 0.0;
	image<rgb>* imgOut = decodeImageWithError(
		K,
		BlockSize,
		encodedBytes.get(),
		outputBytesSize,
		s_quantY, s_quantU, s_quantV,
		dynamicY, dynamicU, dynamicV,
		imgIn, squaredError);
	double mse = squaredError / (imgOut->width() * imgOut->height());
	double psnr = 20.0 * log10(3.0 * 255.0) - 10.0 * log10(mse);
	std::cout << std::format("PSNR {} bpp {} total KB {}", psnr, static_cast<double>(8ULL * outputBytesSize)/static_cast<double>(imgOut->width() * imgOut->height()), outputBytesSize / 1024ULL) << std::endl;
	SaveImageGeneric(imgOut, argv[2], imgFormat::PNG);
	delete imgOut;
	delete imgIn;
	//const auto files = readFiles("D:\\Test\\RAISE");
	//math::Vector patchY(BlockSize * BlockSize);
	//math::Vector patchU(BlockSize * BlockSize);
	//math::Vector patchV(BlockSize * BlockSize);
	//std::vector<BasisChoice> choicesY(K);
	//std::vector<BasisChoice> choicesU(K);
	//std::vector<BasisChoice> choicesV(K);
	//std::mt19937 rand;
	//rand.seed(12345678);
	//const int patchesPerImage = 2000;
	//math::Stat coeffStatsY[K];
	//math::Stat selectStatsY[K];
	//math::Stat coeffStatsU[K];
	//math::Stat selectStatsU[K];
	//math::Stat coeffStatsV[K];
	//math::Stat selectStatsV[K];
	//for (size_t count = 0; count < files.size(); ++count) {
	//	std::string file = files[count];
	//	std::string lowerCaseName(file.size(), 0);
	//	transform(file.cbegin(), file.cend(), lowerCaseName.begin(), ::tolower);
	//	if (lowerCaseName.ends_with(".jpg") || lowerCaseName.ends_with(".tif")) {
	//		std::cout << std::format("processing: {} {}%", file, (100L * count) / files.size()) << std::endl;
	//		imgFormat format;
	//		image<rgb>* image = LoadImageGenericRGB(file.c_str(), &format);
	//		if (image->width() < BlockSize || image->height() < BlockSize) {
	//			delete image;
	//			continue;
	//		}
	//		for (int reps = 0; reps < patchesPerImage; ++reps) {
	//			int x = rand() % (image->width() - BlockSize);
	//			int y = rand() % (image->height() - BlockSize);
	//			for (int offx = 0; offx < BlockSize; ++offx) {
	//				for (int offy = 0; offy < BlockSize; ++offy) {
	//					rgb pt = imRef(image, x + offx, y + offy);
	//					yuv col = YUVFromRGB(pt);
	//					patchY[offx + offy * BlockSize] = col.y;
	//					patchU[offx + offy * BlockSize] = col.u;
	//					patchV[offx + offy * BlockSize] = col.v;
	//				}
	//			}
	//			int countY = CalcMPDynamic(K, s_quantY, choicesY, patchY, dynamicY);
	//			for (int i = 0; i < countY; ++i) {
	//				unsigned short delta = choicesY[i].deltaId;
	//				unsigned short coeff = choicesY[i].intCoeff;
	//				coeffStatsY[i].update((coeff));
	//				selectStatsY[i].update((delta));
	//			}
	//			int countU = CalcMPDynamic(K, s_quantU, choicesU, patchU, dynamicU);
	//			for (int i = 0; i < countU; ++i) {
	//				unsigned short delta = choicesU[i].deltaId;
	//				unsigned short coeff = choicesU[i].intCoeff;
	//				coeffStatsU[i].update((coeff));
	//				selectStatsU[i].update((delta));
	//			}
	//			int countV = CalcMPDynamic(K, s_quantV, choicesV, patchV, dynamicV);
	//			for (int i = 0; i < countV; ++i) {
	//				unsigned short delta = choicesV[i].deltaId;
	//				unsigned short coeff = choicesV[i].intCoeff;
	//				coeffStatsV[i].update((coeff));
	//				selectStatsV[i].update((delta));
	//			}
	//		}
	//		delete image;
	//	}
	//}
	//std::ofstream fileOut;
	//fileOut.open("C:\\Temp\\results6.txt", std::ios::out | std::ios::trunc);
	//fileOut << "Y coeff stats" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("coeff {} min {} max {} mean {} variance {}", i, coeffStatsY[i].min, coeffStatsY[i].max, coeffStatsY[i].mean, coeffStatsY[i].sampleVariance()) << std::endl;
	//}
	//fileOut << "Y basisId stats" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("deltaId {} min {} max {} mean {} variance {}", i, selectStatsY[i].min, selectStatsY[i].max, selectStatsY[i].mean, selectStatsY[i].sampleVariance()) << std::endl;
	//}
	//fileOut << "U coeff stats" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("coeff {} min {} max {} mean {} variance {}", i, coeffStatsU[i].min, coeffStatsU[i].max, coeffStatsU[i].mean, coeffStatsU[i].sampleVariance()) << std::endl;
	//}
	//fileOut << "U basisId stats" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("deltaId {} min {} max {} mean {} variance {}", i, selectStatsU[i].min, selectStatsU[i].max, selectStatsU[i].mean, selectStatsU[i].sampleVariance()) << std::endl;
	//}
	//fileOut << "V coeff stats" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("coeff {} min {} max {} mean {} variance {}", i, coeffStatsV[i].min, coeffStatsV[i].max, coeffStatsV[i].mean, coeffStatsV[i].sampleVariance()) << std::endl;
	//}
	//fileOut << "V basisId stats" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("deltaId {} min {} max {} mean {} variance {}", i, selectStatsV[i].min, selectStatsV[i].max, selectStatsV[i].mean, selectStatsV[i].sampleVariance()) << std::endl;
	//}
	//fileOut.close();
	// 
	//const auto files = readFiles("D:\\Test\\RAISE");
	//math::Vector patchY(BlockSize * BlockSize);
	//math::Vector patchU(BlockSize * BlockSize);
	//math::Vector patchV(BlockSize * BlockSize);
	//std::vector<BasisChoice> choicesY(K);
	//std::vector<BasisChoice> choicesU(K);
	//std::vector<BasisChoice> choicesV(K);
	//std::mt19937 rand;
	//rand.seed(2345678);
	//const int patchesPerImage = 2000;
	//std::vector<std::map<unsigned short, int> > countCY(64);
	//std::vector<std::map<unsigned short, int> > countSY(64);
	//std::vector<std::map<unsigned short, int> > countCU(64);
	//std::vector<std::map<unsigned short, int> > countSU(64);
	//std::vector<std::map<unsigned short, int> > countCV(64);
	//std::vector<std::map<unsigned short, int> > countSV(64);
	//for (size_t count = 0; count < files.size(); ++count) {
	//	std::string file = files[count];
	//	std::string lowerCaseName(file.size(), 0);
	//	transform(file.cbegin(), file.cend(), lowerCaseName.begin(), ::tolower);
	//	if (lowerCaseName.ends_with(".jpg") || lowerCaseName.ends_with(".tif")) {
	//		std::cout << std::format("processing: {} {}%", file, (100L * count) / files.size()) << std::endl;
	//		imgFormat format;
	//		image<rgb>* image = LoadImageGenericRGB(file.c_str(), &format);
	//		if (image->width() < BlockSize || image->height() < BlockSize) {
	//			delete image;
	//			continue;
	//		}
	//		for (int reps = 0; reps < patchesPerImage; ++reps) {
	//			int x = rand() % (image->width() - BlockSize);
	//			int y = rand() % (image->height() - BlockSize);
	//			for (int offx = 0; offx < BlockSize; ++offx) {
	//				for (int offy = 0; offy < BlockSize; ++offy) {
	//					rgb pt = imRef(image, x + offx, y + offy);
	//					yuv col = YUVFromRGB(pt);
	//					patchY[offx + offy * BlockSize] = col.y;
	//					patchU[offx + offy * BlockSize] = col.u;
	//					patchV[offx + offy * BlockSize] = col.v;
	//				}
	//			}
	//			int countY = CalcMPDynamic(K, s_quantY, choicesY, patchY, dynamicY);
	//			for (int i = 0; i < countY; ++i) {
	//				unsigned short delta = choicesY[i].deltaId;
	//				unsigned short coeff = choicesY[i].intCoeff;
	//				countSY[i][delta] = countSY[i][delta] + 1;
	//				countCY[i][coeff] = countCY[i][coeff] + 1;
	//			}
	//			int countU = CalcMPDynamic(K, s_quantU, choicesU, patchU, dynamicU);
	//			for (int i = 0; i < countU; ++i) {
	//				unsigned short delta = choicesU[i].deltaId;
	//				unsigned short coeff = choicesU[i].intCoeff;
	//				countSU[i][delta] = countSU[i][delta] + 1;
	//				countCU[i][coeff] = countCU[i][coeff] + 1;
	//			}
	//			int countV = CalcMPDynamic(K, s_quantV, choicesV, patchV, dynamicV);
	//			for (int i = 0; i < countV; ++i) {
	//				unsigned short delta = choicesV[i].deltaId;
	//				unsigned short coeff = choicesV[i].intCoeff;
	//				countSV[i][delta] = countSV[i][delta] + 1;
	//				countCV[i][coeff] = countCV[i][coeff] + 1;
	//			}
	//		}
	//		delete image;
	//	}
	//}
	//std::vector<int> bestCYM(64);
	//std::vector<int> bestSYM(64);
	//std::vector<int> bestCUM(64);
	//std::vector<int> bestSUM(64);
	//std::vector<int> bestCVM(64);
	//std::vector<int> bestSVM(64);
	//for (int i = 0; i < K; ++i) {
	//	size_t bestCY = 0ULL;
	//	size_t bestSY = 0ULL;
	//	size_t bestCU = 0ULL;
	//	size_t bestSU = 0ULL;
	//	size_t bestCV = 0ULL;
	//	size_t bestSV = 0ULL;
	//	for (int m = 1; m < 4096; ++m) {
	//		size_t bitCountCY = 0ULL;
	//		size_t bitCountSY = 0ULL;
	//		size_t bitCountCU = 0ULL;
	//		size_t bitCountSU = 0ULL;
	//		size_t bitCountCV = 0ULL;
	//		size_t bitCountSV = 0ULL;
	//		for (const auto& freq : countCY[i]) {
	//			bitCountCY += static_cast<size_t>(freq.second * golombCodeLength(freq.first, m));
	//		}
	//		if (m == 1 || bitCountCY < bestCY) {
	//			bestCYM[i] = m;
	//			bestCY = bitCountCY;
	//		}
	//		for (const auto& freq : countSY[i]) {
	//			bitCountSY += static_cast<size_t>(freq.second * golombCodeLength(freq.first, m));
	//		}
	//		if (m == 1 || bitCountSY < bestSY) {
	//			bestSYM[i] = m;
	//			bestSY = bitCountSY;
	//		}
	//		for (const auto& freq : countCU[i]) {
	//			bitCountCU += static_cast<size_t>(freq.second * golombCodeLength(freq.first, m));
	//		}
	//		if (m == 1 || bitCountCU < bestCU) {
	//			bestCUM[i] = m;
	//			bestCU = bitCountCU;
	//		}
	//		for (const auto& freq : countSU[i]) {
	//			bitCountSU += static_cast<size_t>(freq.second * golombCodeLength(freq.first, m));
	//		}
	//		if (m == 1 || bitCountSU < bestSU) {
	//			bestSUM[i] = m;
	//			bestSU = bitCountSU;
	//		}
	//		for (const auto& freq : countCV[i]) {
	//			bitCountCV += static_cast<size_t>(freq.second * golombCodeLength(freq.first, m));
	//		}
	//		if (m == 1 || bitCountCV < bestCV) {
	//			bestCVM[i] = m;
	//			bestCV = bitCountCV;
	//		}
	//		for (const auto& freq : countSV[i]) {
	//			bitCountSV += static_cast<size_t>(freq.second * golombCodeLength(freq.first, m));
	//		}
	//		if (m == 1 || bitCountSV < bestSV) {
	//			bestSVM[i] = m;
	//			bestSV = bitCountSV;
	//		}
	//	}
	//}
	//std::ofstream fileOut;
	//fileOut.open("c:\\temp\\results4.txt", std::ios::out | std::ios::trunc);
	//fileOut <<"static const uint32_t s_golCY[" << K << "] = {" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("{}, ", bestCYM[i]);
	//	if (i % 8 == 7) fileOut << std::endl;
	//}
	//fileOut << "};" << std::endl;
	//fileOut << "static const uint32_t s_golSY[" << K << "] = {" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("{}, ", bestSYM[i]);
	//	if (i % 8 == 7) fileOut << std::endl;
	//}
	//fileOut << "};" << std::endl;

	//fileOut << "static const uint32_t s_golCU[" << K << "] = {" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("{}, ", bestCUM[i]);
	//	if (i % 8 == 7) fileOut << std::endl;
	//}
	//fileOut << "};" << std::endl;

	//fileOut << "static const uint32_t s_golSU[" << K << "] = {" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("{}, ", bestSUM[i]);
	//	if (i % 8 == 7) fileOut << std::endl;
	//}
	//fileOut << "};" << std::endl;

	//fileOut << "static const uint32_t s_golCV[" << K << "] = {" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("{}, ", bestCVM[i]);
	//	if (i % 8 == 7) fileOut << std::endl;
	//}
	//fileOut << "};" << std::endl;

	//fileOut << "static const uint32_t s_golSV[" << K << "] = {" << std::endl;
	//for (int i = 0; i < K; ++i) {
	//	fileOut << std::format("{}, ", bestSVM[i]);
	//	if (i % 8 == 7) fileOut << std::endl;
	//}
	//fileOut << "};" << std::endl;

	//fileOut.close();
	// 
	//const auto files = readFiles("D:\\Test\\RAISE");
	//math::Vector patchY(BlockSize* BlockSize);
	//std::vector<BasisChoice> choicesY(K);
	//std::mt19937 rand;
	//rand.seed(2345678);
	//std::vector<double> shift(31);
	//std::vector<double> peak(31);
	//for (int i = 0; i < 31; ++i) {
	//	shift[i] = 0.5;
	//}
	//for (int iter = 0;iter < 1000;++iter) {
	//	for (int choice = 1; choice < 31; ++choice) {
	//		while (true) {
	//			double totalError = 0.0;
	//			double N = 0.0;
	//			s_quantY[choice] = s_quantY[choice] + shift[choice - 1];
	//			for (int samples = 0; samples < 10; ++samples) {
	//				std::string file = files[rand() % files.size()];
	//				std::string lowerCaseName(file.size(), 0);
	//				transform(file.cbegin(), file.cend(), lowerCaseName.begin(), ::tolower);
	//				if (!(lowerCaseName.ends_with(".jpg") || lowerCaseName.ends_with(".tif"))) {
	//					continue;
	//				}
	//				imgFormat format;
	//				image<double>* image = LoadImageGenericMono(file.c_str(), &format);
	//				if (image->width() < BlockSize || image->height() < BlockSize) {
	//					delete image;
	//					continue;
	//				}
	//				for (int reps = 0; reps < 200; ++reps) {
	//					int x = rand() % (image->width() - BlockSize);
	//					int y = rand() % (image->height() - BlockSize);
	//					for (int offx = 0; offx < BlockSize; ++offx) {
	//						for (int offy = 0; offy < BlockSize; ++offy) {
	//							double pt = imRef(image, x + offx, y + offy);
	//							patchY[offx + offy * BlockSize] = pt;
	//						}
	//					}
	//					int countY = CalcMPDynamic(K, s_quantY, choicesY, patchY, dynamicY);
	//					math::Vector decodedY = FromCoeffsDynamic(countY, s_quantY, choicesY, dynamicY);
	//					for (int offx = 0; offx < BlockSize; ++offx) {
	//						for (int offy = 0; offy < BlockSize; ++offy) {
	//							double pt = imRef(image, x + offx, y + offy);
	//							double err = decodedY[offx + offy * BlockSize] - pt;
	//							totalError += square(err);
	//							N += 1.0;
	//						}
	//					}
	//				}
	//				delete image;
	//			}
	//			double mse = totalError / N;
	//			double psnr = 20.0 * log10(255.0) - 10.0 * log10(mse);
	//			if (psnr < 50.0) {
	//				std::cout << std::format("fail choice {} at {} with psnr {}", choice, s_quantY[choice], psnr) << std::endl;
	//				s_quantY[choice] = std::max(s_quantY[choice] - shift[choice - 1], 1.0);
	//				shift[choice - 1] = std::max(shift[choice - 1] / 2.0, 0.125);
	//				break;
	//			} else {
	//				for (int i = 0; i < K; ++i) {
	//					std::cout << s_quantY[i] << ", ";
	//				}
	//				std::cout << psnr << std::endl;
	//			}
	//		}
	//	}
	//}
	return 0;
}

