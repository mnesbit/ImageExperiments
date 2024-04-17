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
#include "MatchingPursuit.h"
#include "BitBuffer.h"
#include "BasisSet.h"
#include "../ImageHelper/inc/imgloader.h"
#include "../SimpleMatrix/inc/covariance.h"
#include "../SimpleMatrix/inc/symmeigen.h"
#include "../ImageHelper/inc/misc.h"

using namespace img;
using namespace matching;
using namespace bitbuffer;
using namespace basis;

static double s_quantY[32] = {
	4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
};

static const double s_quantU[32] = {
	4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
};

static const double s_quantV[32] = {
	4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
};

static const uint32_t s_golCY[32] = {
511, 63, 31, 31, 15, 15, 15, 15,
15, 15, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7,
3, 3, 3, 3, 3, 3, 3, 3,
};

static const uint32_t s_golSY[32] = {
1, 511, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
};

static const uint32_t s_golCU[32] = {
63, 15, 7, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3,
};

static const uint32_t s_golSU[32] = {
1, 511, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
};

static const uint32_t s_golCV[32] = {
31, 15, 7, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3,
};

static const uint32_t s_golSV[32] = {
7, 511, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255,
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

struct Stat {
	double N;
	double min;
	double max;
	double mean;
	double sumSq;

	void update(double val) {
		if (N == 0) {
			N = 1.0;
			min = val;
			max = val;
			mean = val;
			sumSq = 0.0;
		} else {
			if (val < min) min = val;
			if (val > max) max = val;
			N += 1.0;
			double delta = val - mean;
			mean += delta / N;
			double delta2 = val - mean;
			sumSq += delta * delta2;
		}
	}

	double sampleVariance() const {
		return sumSq / (N - 1.0);
	}
};

std::unique_ptr<uint8_t[]> encodeImage(const image<rgb>* imgIn, const size_t K, const size_t blockSize,
	DynamicDictionaryFunction dynamicY, DynamicDictionaryFunction dynamicU, DynamicDictionaryFunction dynamicV,
	size_t& outputByteSize) {
	const size_t width = imgIn->width();
	const size_t height = imgIn->height();
	BitBuffer bitsOut;
	bitsOut.WriteInt(static_cast<uint32_t>(width));
	bitsOut.WriteInt(static_cast<uint32_t>(height));
	std::vector<BasisChoice> choicesY(K);
	std::vector<BasisChoice> choicesU(K);
	std::vector<BasisChoice> choicesV(K);
	math::Vector blockY(blockSize * blockSize);
	math::Vector blockU(blockSize * blockSize);
	math::Vector blockV(blockSize * blockSize);
	for (size_t x = 0; x < width; x += blockSize) {
		std::cout << std::format("encode {} %", (100 * x) / width) << std::endl;
		for (size_t y = 0; y < height; y += blockSize) {
			for (size_t dx = 0; dx < blockSize; ++dx) {
				size_t u = x + dx;
				for (size_t dy = 0; dy < blockSize; ++dy) {
					size_t v = y + dy;
					if ((u < width) && (v < height)) {
						rgb pt = imRef(imgIn, u, v);
						yuv col = YUVFromRGB(pt);
						blockY[dx + (blockSize * dy)] = col.y;
						blockU[dx + (blockSize * dy)] = col.u;
						blockV[dx + (blockSize * dy)] = col.v;
					} else {
						blockY[dx + (blockSize * dy)] = 0.0;
						blockU[dx + (blockSize * dy)] = 0.0;
						blockV[dx + (blockSize * dy)] = 0.0;
					}
				}
			}
			size_t countY = CalcMPDynamic(K, s_quantY, choicesY, blockY, dynamicY);
			bitsOut.WriteBits(countY, std::bit_width(K));
			if (countY > 0) {
				if (choicesY[0].deltaId == 0) { // almost always DC components first
					bitsOut.WriteBits(0, 1);
					writeGolombCode(choicesY[0].intCoeff, s_golCY[0], bitsOut);
				} else {
					bitsOut.WriteBits(1, 1);
					writeGolombCode(choicesY[0].deltaId, s_golSY[0], bitsOut);
					writeGolombCode(choicesY[0].intCoeff, s_golCY[0], bitsOut);
				}
				for (int i = 1; i < countY; ++i) {
					writeGolombCode(choicesY[i].deltaId, s_golSY[i], bitsOut);
					writeGolombCode(choicesY[i].intCoeff, s_golCY[i], bitsOut);
				}
			}
			size_t countU = CalcMPDynamic(K, s_quantU, choicesU, blockU, dynamicU);
			bitsOut.WriteBits(countU, std::bit_width(K));
			if (countU > 0) {
				if (choicesU[0].deltaId == 0) { // almost always DC components first
					bitsOut.WriteBits(0, 1);
					writeGolombCode(choicesU[0].intCoeff, s_golCU[0], bitsOut);
				} else {
					bitsOut.WriteBits(1, 1);
					writeGolombCode(choicesU[0].deltaId, s_golSU[0], bitsOut);
					writeGolombCode(choicesU[0].intCoeff, s_golCU[0], bitsOut);
				}
				for (int i = 1; i < countU; ++i) {
					writeGolombCode(choicesU[i].deltaId, s_golSU[i], bitsOut);
					writeGolombCode(choicesU[i].intCoeff, s_golCU[i], bitsOut);
				}
			}
			size_t countV = CalcMPDynamic(K, s_quantV, choicesV, blockV, dynamicV);
			bitsOut.WriteBits(countV, std::bit_width(K));
			if (countV > 0) {
				if (choicesV[0].deltaId == 0) { // almost always DC components first
					bitsOut.WriteBits(0, 1);
					writeGolombCode(choicesV[0].intCoeff, s_golCV[0], bitsOut);
				} else {
					bitsOut.WriteBits(1, 1);
					writeGolombCode(choicesV[0].deltaId, s_golSV[0], bitsOut);
					writeGolombCode(choicesV[0].intCoeff, s_golCV[0], bitsOut);
				}
				for (int i = 1; i < countV; ++i) {
					writeGolombCode(choicesV[i].deltaId, s_golSV[i], bitsOut);
					writeGolombCode(choicesV[i].intCoeff, s_golCV[i], bitsOut);
				}
			}
		}
	}
	return bitsOut.Save(outputByteSize);
}

image<rgb>* decodeImage(const size_t K, const size_t blockSize,const uint8_t bytes[], size_t byteSize,
	DynamicDictionaryFunction dynamicY, DynamicDictionaryFunction dynamicU, DynamicDictionaryFunction dynamicV,
	const image<rgb>* original, double& squaredError) {
	BitBuffer bitsIn;
	size_t lengthsCY[64];
	size_t lengthsCU[64];
	size_t lengthsCV[64];
	size_t lengthsSY[64];
	size_t lengthsSU[64];
	size_t lengthsSV[64];
	for (int i = 0; i < 64; ++i) {
		lengthsCY[i] = 0;
		lengthsCU[i] = 0;
		lengthsCV[i] = 0;
		lengthsSY[i] = 0;
		lengthsSU[i] = 0;
		lengthsSV[i] = 0;
	}
	bitsIn.Load(bytes, 0ULL, 8ULL * byteSize);
	const size_t width = static_cast<size_t>(bitsIn.ReadInt());
	const size_t height = static_cast<size_t>(bitsIn.ReadInt());
	image<rgb>* imgOut = new image<rgb>(width, height, false);
	std::vector<BasisChoice> choicesY(K);
	std::vector<BasisChoice> choicesU(K);
	std::vector<BasisChoice> choicesV(K);
	for (size_t x = 0; x < width; x += blockSize) {
		std::cout << std::format("decode {} %", (100 * x) / width) << std::endl;
		for (size_t y = 0; y < height; y += blockSize) {
			size_t countY = bitsIn.ReadBits(std::bit_width(K));
			lengthsSY[0] += std::bit_width(K);
			if (countY > 0) {
				if (bitsIn.ReadBits(1) == 0) {
					choicesY[0].deltaId = 0;
					choicesY[0].intCoeff = readGolombCode(s_golCY[0], bitsIn);
					lengthsSY[0] += 1;
					lengthsCY[0] += golombCodeLength(choicesY[0].intCoeff, s_golCY[0]);
				} else {
					choicesY[0].deltaId = readGolombCode(s_golSY[0], bitsIn);
					choicesY[0].intCoeff = readGolombCode(s_golCY[0], bitsIn);
					lengthsSY[0] += golombCodeLength(choicesY[0].deltaId, s_golSY[0]) + 1;
					lengthsCY[0] += golombCodeLength(choicesY[0].intCoeff, s_golCY[0]);
				}
				for (int i = 1; i < countY; ++i) {
					choicesY[i].deltaId = readGolombCode(s_golSY[i], bitsIn);
					choicesY[i].intCoeff = readGolombCode(s_golCY[i], bitsIn);
					lengthsSY[i] += golombCodeLength(choicesY[i].deltaId, s_golSY[i]);
					lengthsCY[i] += golombCodeLength(choicesY[i].intCoeff, s_golCY[i]);
				}
			}
			math::Vector decodedY = FromCoeffsDynamic(countY, s_quantY, choicesY, dynamicY);
			size_t countU = bitsIn.ReadBits(std::bit_width(K));
			lengthsSU[0] += std::bit_width(K);
			if (countU > 0) {
				if (bitsIn.ReadBits(1) == 0) {
					choicesU[0].deltaId = 0;
					choicesU[0].intCoeff = readGolombCode(s_golCU[0], bitsIn);
					lengthsSU[0] += 1;
					lengthsCU[0] += golombCodeLength(choicesU[0].intCoeff, s_golCU[0]);
				} else {
					choicesU[0].deltaId = readGolombCode(s_golSU[0], bitsIn);
					choicesU[0].intCoeff = readGolombCode(s_golCU[0], bitsIn);
					lengthsSU[0] += golombCodeLength(choicesU[0].deltaId, s_golSU[0]) + 1;
					lengthsCU[0] += golombCodeLength(choicesU[0].intCoeff, s_golCU[0]);
				}
				for (int i = 1; i < countU; ++i) {
					choicesU[i].deltaId = readGolombCode(s_golSU[i], bitsIn);
					choicesU[i].intCoeff = readGolombCode(s_golCU[i], bitsIn);
					lengthsSU[i] += golombCodeLength(choicesU[i].deltaId, s_golSU[i]);
					lengthsCU[i] += golombCodeLength(choicesU[i].intCoeff, s_golCU[i]);
				}
			}
			math::Vector decodedU = FromCoeffsDynamic(countU, s_quantU, choicesU, dynamicU);
			size_t countV = bitsIn.ReadBits(std::bit_width(K));
			lengthsSV[0] += std::bit_width(K);
			if (countV > 0) {
				if (bitsIn.ReadBits(1) == 0) {
					choicesV[0].deltaId = 0;
					choicesV[0].intCoeff = readGolombCode(s_golCV[0], bitsIn);
					lengthsSV[0] += 1;
					lengthsCV[0] += golombCodeLength(choicesV[0].intCoeff, s_golCV[0]);
				} else {
					choicesV[0].deltaId = readGolombCode(s_golSV[0], bitsIn);
					choicesV[0].intCoeff = readGolombCode(s_golCV[0], bitsIn);
					lengthsSV[0] += golombCodeLength(choicesV[0].deltaId, s_golSV[0]) + 1;
					lengthsCV[0] += golombCodeLength(choicesV[0].intCoeff, s_golCV[0]);
				}
				for (int i = 1; i < countV; ++i) {
					choicesV[i].deltaId = readGolombCode(s_golSV[i], bitsIn);
					choicesV[i].intCoeff = readGolombCode(s_golCV[i], bitsIn);
					lengthsSV[i] += golombCodeLength(choicesV[i].deltaId, s_golSV[i]);
					lengthsCV[i] += golombCodeLength(choicesV[i].intCoeff, s_golCV[i]);
				}
			}
			math::Vector decodedV = FromCoeffsDynamic(countV, s_quantV, choicesV, dynamicV);
			for (size_t dx = 0; dx < blockSize; ++dx) {
				size_t u = x + dx;
				for (size_t dy = 0; dy < blockSize; ++dy) {
					size_t v = y + dy;
					if ((u < width) && (v < height)) {
						rgb ptOrig = imRef(original, u, v);
						rgb pt = RGBFromYUV(yuv{
							.y = decodedY[dx + (blockSize * dy)],
							.u = decodedU[dx + (blockSize * dy)],
							.v = decodedV[dx + (blockSize * dy)]
							});
						imRef(imgOut, u, v) = pt;
						double errR = static_cast<double>(ptOrig.r) - static_cast<double>(pt.r);
						double errG = static_cast<double>(ptOrig.g) - static_cast<double>(pt.g);
						double errB = static_cast<double>(ptOrig.b) - static_cast<double>(pt.b);
						squaredError += square(errR) + square(errG) + square(errB);
					}
				}
			}
		}
	}
	for (int i = 0; i < 64; ++i) {
		std::cout << std::format("index {} Y bits {}/{} U bits {}/{} V bits {}/{}", i, lengthsSY[i], lengthsCY[i], lengthsSU[i], lengthsCU[i], lengthsSV[i], lengthsCV[i]) << std::endl;
	}
	return imgOut;
}

int main(int argc, char* argv[]) {
	if (argc != 3) {
		std::cout << "usage: Compression.exe <input file> <output file>" << std::endl;
		return -1;
	}
	const int BlockSize = 8;
	const int K = 32;
	std::cout << "Calculating basis vectors" << std::endl;
	std::vector<std::vector<double> > basisSet = distinctLineShapes(BlockSize);
	math::Matrix baseDict = createSegmentDictionary(BlockSize, basisSet);
	size_t fullRowCount = baseDict.Rows();
	std::vector<math::Matrix> detailBasisY;
	std::vector<math::Matrix> detailBasisU;
	std::vector<math::Matrix> detailBasisV;
	for (const std::vector<double>& basis : basisSet) {
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
		int offset = baseDict.Rows();
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
				offset += static_cast<size_t>(detail.Rows());
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
	std::unique_ptr<uint8_t[]> encodedBytes = encodeImage(imgIn, K, BlockSize, dynamicY, dynamicU, dynamicV, outputBytesSize);
	double squaredError = 0.0;
	image<rgb>* imgOut = decodeImage(K, BlockSize, encodedBytes.get(), outputBytesSize, dynamicY, dynamicU, dynamicV, imgIn, squaredError);
	double mse = squaredError / (imgOut->width() * imgOut->height());
	double psnr = 20.0 * log10(3.0 * 255.0) - 10.0 * log10(mse);
	std::cout << std::format("PSNR {} bpp {} total KB {}", psnr, static_cast<double>(8ULL * outputBytesSize)/static_cast<double>(imgOut->width() * imgOut->height()), outputBytesSize / 1024ULL) << std::endl;
	SaveImageGeneric(imgOut, argv[2], imgFormat::PNG);
	delete imgOut;
	//const auto files = readFiles(argv[1]);
	//math::Vector patchY(BlockSize * BlockSize);
	//math::Vector patchU(BlockSize * BlockSize);
	//math::Vector patchV(BlockSize * BlockSize);
	//std::vector<BasisChoice> choicesY(K);
	//std::vector<BasisChoice> choicesU(K);
	//std::vector<BasisChoice> choicesV(K);
	//std::mt19937 rand;
	//rand.seed(2345678);
	//const int patchesPerImage = 2000;
	//Stat coeffStatsY[K];
	//Stat selectStatsY[K];
	//Stat coeffStatsU[K];
	//Stat selectStatsU[K];
	//Stat coeffStatsV[K];
	//Stat selectStatsV[K];
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
	//			CalcMPDynamic(K, s_quantY, choicesY, patchY, dynamicY);
	//			for (int i = 0; i < K; ++i) {
	//				unsigned short delta = choicesY[i].deltaId;
	//				unsigned short coeff = choicesY[i].intCoeff;
	//				coeffStatsY[i].update((coeff));
	//				selectStatsY[i].update((delta));
	//				if (delta == 0 && coeff == 0) {
	//					break;
	//				}
	//			}
	//			CalcMPDynamic(K, s_quantU, choicesU, patchU, dynamicU);
	//			for (int i = 0; i < K; ++i) {
	//				unsigned short delta = choicesU[i].deltaId;
	//				unsigned short coeff = choicesU[i].intCoeff;
	//				coeffStatsU[i].update((coeff));
	//				selectStatsU[i].update((delta));
	//				if (delta == 0 && coeff == 0) {
	//					break;
	//				}
	//			}
	//			CalcMPDynamic(K, s_quantV, choicesV, patchV, dynamicV);
	//			for (int i = 0; i < K; ++i) {
	//				unsigned short delta = choicesV[i].deltaId;
	//				unsigned short coeff = choicesV[i].intCoeff;
	//				coeffStatsV[i].update((coeff));
	//				selectStatsV[i].update((delta));
	//				if (delta == 0 && coeff == 0) {
	//					break;
	//				}
	//			}
	//		}
	//		delete image;
	//	}
	//}
	//std::ofstream fileOut;
	//fileOut.open(argv[2], std::ios::out|std::ios::trunc);
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
	//std::vector<double> shift(63);
	//std::vector<double> peak(63);
	//for (int i = 0; i < 63; ++i) {
	//	shift[i] = 0.5;
	//}
	//for (int iter = 0;iter < 1000;++iter) {
	//	for (int choice = 1; choice < 64; ++choice) {
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
	//					CalcMPDynamic(K, s_quantY, choicesY, patchY, dynamicY);
	//					math::Vector decodedY = FromCoeffsDynamic(K, s_quantY, choicesY, dynamicY);
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
	//			if (psnr < 45.0) {
	//				std::cout << std::format("fail choice {} at {} with psnr {}", choice, s_quantY[choice], psnr) << std::endl;
	//				s_quantY[choice] = std::max(s_quantY[choice] - shift[choice - 1], 1.0);
	//				shift[choice - 1] = std::max(shift[choice - 1] / 2.0, 0.125);
	//				break;
	//			} else {
	//				for (int i = 0; i < BlockSize * BlockSize; ++i) {
	//					std::cout << s_quantY[i] << ", ";
	//				}
	//				std::cout << psnr << std::endl;
	//			}
	//		}
	//	}
	//}
	return 0;
}

