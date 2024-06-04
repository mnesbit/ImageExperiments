#include <iostream>
#include <cstdio>
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

const int BlockSize = 8;
const int K = 32;

static std::vector<std::string> readFiles(std::string path) {
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

static void writeDictionaryToPNG(const math::Matrix& dictionary, std::string fileName) {
	size_t blockSize = static_cast<size_t>(sqrt(dictionary.Columns()));
	if (dictionary.Columns() != blockSize * blockSize) {
		throw std::invalid_argument("Dictionary should have square Column count");
	}
	size_t approxWidth = static_cast<size_t>(sqrt(dictionary.Rows()));
	size_t approxHeight = (dictionary.Rows() / approxWidth) + ((dictionary.Rows() % approxWidth == 0) ? 0 : 1);
	std::unique_ptr<image<double> > basisPic = std::make_unique<image<double> >(approxWidth * blockSize, approxHeight * blockSize);
	for (size_t i = 0; i < dictionary.Rows(); ++i) {
		size_t blockx = i % approxWidth;
		size_t blocky = i / approxWidth;
		for (size_t x = 0; x < blockSize; ++x) {
			for (size_t y = 0; y < blockSize; ++y) {
				imRef(basisPic, x + blockx * blockSize, y + blocky * blockSize) = dictionary[i][x + (y * blockSize)];
			}
		}
	}
	SaveImageGeneric(basisPic.get(), fileName.c_str(), imgFormat::PNG);
	basisPic.release();
}

static void printUsage() {
	std::cout << "To compress to file" << std::endl;
	std::cout << "usage: Compression.exe -c <quality 1.0-8.0> <input file path> <output compressed file path>" << std::endl;
	std::cout << "To decompress to file" << std::endl;
	std::cout << "usage: Compression.exe -d <input compressed file path> <output PNG file path>" << std::endl;
	std::cout << "To display psnr and bpp of compression cycle" << std::endl;
	std::cout << "usage: Compression.exe -n  <quality 1.0-8.0> <input compressed file path> [optional output PNG path]" << std::endl;
	std::cout << "To display psnr and bpp of JPEG compression cycle at quality R" << std::endl;
	std::cout << "usage: Compression.exe -j <quality 1-100> <input compressed file path>" << std::endl;
	std::cout << "Enumerate all .TIF/.JPG files in folder and gather statistics (to use as basis of future quantization calcs)" << std::endl;
	std::cout << "usage: Compression.exe -s <seed> <patches per image> <input folder path> <outputs statistics file path>" << std::endl;
	std::cout << "Enumerate all .TIF/.JPG files in folder and calclulate compression curves for them" << std::endl;
	std::cout << "Compression.exe -g <input folder path> <outputs statistics file path>" << std::endl;
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		printUsage();
		return -1;
	}
	std::string mode = argv[1];
	if (mode == "-c") {
		if (argc != 5) {
			printUsage();
			return -1;
		}
		std::cout << "create basis function" << std::endl;
		std::unique_ptr<image<rgb> > imgIn = LoadImageGenericRGB(argv[3]);
		std::unique_ptr<CompressionContextFast> context;
		if (std::string(argv[2]) == "max") {
			context = createCompressionContextFast(K, BlockSize, 0.0);
			for (size_t i = 0; i < K; ++i) {
				context->Y.Quant[i] = 1.0;
				context->U.Quant[i] = 1.0;
				context->V.Quant[i] = 1.0;
			}
		} else {
			double bppAllocation = atof(argv[2]);
			context = createCompressionContextFast(K, BlockSize, bppAllocation);
		}
		std::cout << std::format("start processing image {} using largest {} coefficients.", argv[3], K) << std::endl;
		size_t outputBytesSize;
		std::unique_ptr<uint8_t[]> encodedBytes = encodeImageFast(
			imgIn.get(),
			K,
			BlockSize,
			context->Y.Quant, context->U.Quant, context->V.Quant,
			context->Y.Dynamic, context->U.Dynamic, context->V.Dynamic,
			outputBytesSize);
		std::cout << std::format("bpp {} total KB {}", static_cast<double>(8ULL * outputBytesSize) / static_cast<double>(imgIn->width() * imgIn->height()), outputBytesSize / 1024ULL) << std::endl;
		imgIn.reset();
		std::ofstream fileOut;
		fileOut.open(argv[4], std::ios::out | std::ios::trunc |std::ios::binary);
		fileOut.write(reinterpret_cast<const char*>(encodedBytes.get()), outputBytesSize);
		fileOut.close();
	} else if (mode == "-d") {
		if (argc != 4) {
			printUsage();
			return -1;
		}
		std::ifstream fileIn(argv[2], std::ios::binary |std::ios::in | std::ios::ate);
		std::ifstream::pos_type len = fileIn.tellg();
		std::unique_ptr<uint8_t[]> buffer = std::make_unique<uint8_t[]>(len);
		fileIn.seekg(0, std::ios::beg);
		fileIn.read(reinterpret_cast<char*>(buffer.get()), len);
		std::unique_ptr<image<rgb> > imgOut = decodeImageFast(buffer.get(), len);
		buffer.reset();
		img::SaveImageGeneric(imgOut.get(), argv[3], imgFormat::PNG);
		imgOut.reset();
	} else if (mode == "-n") {
		if (argc != 4 && argc != 5) {
			printUsage();
			return -1;
		}
		const auto start = std::chrono::system_clock::now();
		std::cout << "create basis function" << std::endl;
		std::unique_ptr<image<rgb> > imgIn = LoadImageGenericRGB(argv[3]);
		std::unique_ptr<CompressionContextFast> context;
		if (std::string(argv[2]) == "max") {
			context = createCompressionContextFast(K, BlockSize, 0.0);
			for (size_t i = 0; i < K; ++i) {
				context->Y.Quant[i] = 1.0;
				context->U.Quant[i] = 1.0;
				context->V.Quant[i] = 1.0;
			}
		} else {
			double bppAllocation = atof(argv[2]);
			context = createCompressionContextFast(K, BlockSize, bppAllocation);
		}
		std::cout << std::format("start processing image {} using largest {} coefficients.", argv[3], K) << std::endl;
		size_t outputBytesSize;
		std::unique_ptr<uint8_t[]> encodedBytes = encodeImageFast(
			imgIn.get(),
			K,
			BlockSize,
			context->Y.Quant, context->U.Quant, context->V.Quant,
			context->Y.Dynamic, context->U.Dynamic, context->V.Dynamic,
			outputBytesSize);
		std::unique_ptr<image<rgb> > imgOut = decodeImageFast(encodedBytes.get(), outputBytesSize);
		const auto end = std::chrono::system_clock::now();
		std::cout << duration_cast<std::chrono::milliseconds>(end - start) << std::endl;
		double psnr = calculatePSNR(imgIn.get(), imgOut.get());
		std::cout << std::format("PSNR {} bpp {} total KB {}", psnr, static_cast<double>(8ULL * outputBytesSize) / static_cast<double>(imgOut->width() * imgOut->height()), outputBytesSize / 1024ULL) << std::endl;
		imgIn.reset();
		if (argc == 5) {
			SaveImageGeneric(imgOut.get(), argv[4], imgFormat::PNG);
		}
		imgOut.reset();
	} else if (mode == "-j") {
		if (argc != 4) {
			printUsage();
			return -1;
		}
		int quality = atoi(argv[2]);
		std::unique_ptr<image<rgb> > imgIn = LoadImageGenericRGB(argv[3]);
		char tmpFileName[L_tmpnam_s];
		tmpnam_s(tmpFileName, L_tmpnam_s);
		img::SaveImageGeneric(imgIn.get(), tmpFileName, imgFormat::JPEG, quality); //convert to JPEG
		size_t jpegSize = std::filesystem::file_size(tmpFileName);
		std::unique_ptr<image<rgb> > imgConv = LoadImageGenericRGB(tmpFileName);
		double jpegpsnr = calculatePSNR(imgIn.get(), imgConv.get());
		imgConv.reset();
		std::cout << std::format("JPEG PSNR {} bpp {} total KB {}", jpegpsnr, static_cast<double>(8ULL * jpegSize) / static_cast<double>(imgIn->width() * imgIn->height()), jpegSize / 1024ULL) << std::endl;
		imgIn.reset();
		std::filesystem::remove(tmpFileName);
	} else if (mode == "-s") {
		if (argc != 6) {
			printUsage();
		}
		const auto files = readFiles(argv[4]);
		math::Vector patchY(BlockSize * BlockSize);
		math::Vector patchU(BlockSize * BlockSize);
		math::Vector patchV(BlockSize * BlockSize);
		std::vector<BasisChoice> choicesY(K);
		std::vector<BasisChoice> choicesU(K);
		std::vector<BasisChoice> choicesV(K);
		std::mt19937 rand;
		rand.seed(atoi(argv[2]));
		const int patchesPerImage = atoi(argv[3]);
		math::Stat coeffStatsY[K]{};
		math::Stat selectStatsY[K]{};
		math::Stat coeffStatsU[K]{};
		math::Stat selectStatsU[K]{};
		math::Stat coeffStatsV[K]{};
		math::Stat selectStatsV[K]{};
		std::unique_ptr<CompressionContext> context = createCompressionContext(K, BlockSize, 0.0);
		for (size_t i = 0; i < K; ++i) {
			context->Y.Quant[i] = 1.0;
			context->U.Quant[i] = 1.0;
			context->V.Quant[i] = 1.0;
		}
		for (size_t count = 0; count < files.size(); ++count) {
			std::string file = files[count];
			std::string lowerCaseName(file.size(), 0);
			transform(file.cbegin(), file.cend(), lowerCaseName.begin(), ::tolower);
			if (lowerCaseName.ends_with(".jpg") || lowerCaseName.ends_with(".tif")) {
				std::cout << std::format("processing: {} {}%", file, (100L * count) / files.size()) << std::endl;
				imgFormat format;
				std::unique_ptr<image<rgb> > image = LoadImageGenericRGB(file.c_str(), &format);
				if (image->width() < BlockSize || image->height() < BlockSize) {
					image.reset();
					continue;
				}
				for (int reps = 0; reps < patchesPerImage; ++reps) {
					int x = rand() % (image->width() - BlockSize);
					int y = rand() % (image->height() - BlockSize);
					for (int offx = 0; offx < BlockSize; ++offx) {
						for (int offy = 0; offy < BlockSize; ++offy) {
							rgb pt = imRef(image, x + offx, y + offy);
							yuv col = YUVFromRGB(pt);
							patchY[offx + offy * BlockSize] = col.y;
							patchU[offx + offy * BlockSize] = col.u;
							patchV[offx + offy * BlockSize] = col.v;
						}
					}
					int countY = CalcMPDynamic(K, context->Y.Quant.Data(), choicesY, patchY, context->Y.Dynamic);
					for (int i = 0; i < countY; ++i) {
						unsigned short delta = choicesY[i].deltaId;
						unsigned short coeff = choicesY[i].intCoeff;
						coeffStatsY[i].update((coeff));
						selectStatsY[i].update((delta));
					}
					int countU = CalcMPDynamic(K, context->U.Quant.Data(), choicesU, patchU, context->U.Dynamic);
					for (int i = 0; i < countU; ++i) {
						unsigned short delta = choicesU[i].deltaId;
						unsigned short coeff = choicesU[i].intCoeff;
						coeffStatsU[i].update((coeff));
						selectStatsU[i].update((delta));
					}
					int countV = CalcMPDynamic(K, context->V.Quant.Data(), choicesV, patchV, context->V.Dynamic);
					for (int i = 0; i < countV; ++i) {
						unsigned short delta = choicesV[i].deltaId;
						unsigned short coeff = choicesV[i].intCoeff;
						coeffStatsV[i].update((coeff));
						selectStatsV[i].update((delta));
					}
				}
				image.reset();
			}
		}
		std::ofstream fileOut;
		fileOut.open(argv[5], std::ios::out | std::ios::trunc);
		fileOut << "Y coeff stats" << std::endl;
		for (int i = 0; i < K; ++i) {
			fileOut << std::format("coeff {} min {} max {} range {} mean {} variance {} std dev {}", i, coeffStatsY[i].min, coeffStatsY[i].max, coeffStatsY[i].max - coeffStatsY[i].min, coeffStatsY[i].mean, coeffStatsY[i].sampleVariance(), sqrt(coeffStatsY[i].sampleVariance())) << std::endl;
		}
		fileOut << "Y basisId stats" << std::endl;
		for (int i = 0; i < K; ++i) {
			fileOut << std::format("deltaId {} min {} max {} range {} mean {} variance {} std dev {}", i, selectStatsY[i].min, selectStatsY[i].max, selectStatsY[i].max - selectStatsY[i].min, selectStatsY[i].mean, selectStatsY[i].sampleVariance(), sqrt(selectStatsY[i].sampleVariance())) << std::endl;
		}
		fileOut << "U coeff stats" << std::endl;
		for (int i = 0; i < K; ++i) {
			fileOut << std::format("coeff {} min {} max {} range {} mean {} variance {} std dev {}", i, coeffStatsU[i].min, coeffStatsU[i].max, coeffStatsU[i].max - coeffStatsU[i].min, coeffStatsU[i].mean, coeffStatsU[i].sampleVariance(), sqrt(coeffStatsU[i].sampleVariance())) << std::endl;
		}
		fileOut << "U basisId stats" << std::endl;
		for (int i = 0; i < K; ++i) {
			fileOut << std::format("deltaId {} min {} max {} range {} mean {} variance {} std dev {}", i, selectStatsU[i].min, selectStatsU[i].max, selectStatsU[i].max - selectStatsU[i].min, selectStatsU[i].mean, selectStatsU[i].sampleVariance(), sqrt(selectStatsU[i].sampleVariance())) << std::endl;
		}
		fileOut << "V coeff stats" << std::endl;
		for (int i = 0; i < K; ++i) {
			fileOut << std::format("coeff {} min {} max {} range {} mean {} variance {} std dev {}", i, coeffStatsV[i].min, coeffStatsV[i].max, coeffStatsV[i].max - coeffStatsV[i].min, coeffStatsV[i].mean, coeffStatsV[i].sampleVariance(), sqrt(coeffStatsV[i].sampleVariance())) << std::endl;
		}
		fileOut << "V basisId stats" << std::endl;
		for (int i = 0; i < K; ++i) {
			fileOut << std::format("deltaId {} min {} max {} range {} mean {} variance {} std dev {}", i, selectStatsV[i].min, selectStatsV[i].max, selectStatsV[i].max - selectStatsV[i].min, selectStatsV[i].mean, selectStatsV[i].sampleVariance(), sqrt(selectStatsV[i].sampleVariance())) << std::endl;
		}
		fileOut.close();
		return -1;
	} else if (mode == "-g") {
		if (argc != 4) {
			printUsage();
		}
		const auto files = readFiles(argv[2]);
		std::ofstream fileOut;
		fileOut.open(argv[3], std::ios::out | std::ios::trunc);
		fileOut << "File, Mode, Quality, Size, BPP, PSNR" << std::endl;
		for (size_t count = 0; count < files.size(); ++count) {
			std::string file = files[count];
			std::string lowerCaseName(file.size(), 0);
			transform(file.cbegin(), file.cend(), lowerCaseName.begin(), ::tolower);
			if (lowerCaseName.ends_with(".jpg") || lowerCaseName.ends_with(".tif")) {
				std::cout << std::format("JPEG processing: {} {}%", file, (100L * count) / files.size()) << std::endl;
				std::unique_ptr<image<rgb> > imgIn = LoadImageGenericRGB(file.c_str());
				for (int quality = 100; quality >= 10; quality -= 20) {
					char tmpFileName[L_tmpnam_s];
					tmpnam_s(tmpFileName, L_tmpnam_s);
					img::SaveImageGeneric(imgIn.get(), tmpFileName, imgFormat::JPEG, quality); //convert to JPEG
					size_t jpegSize = std::filesystem::file_size(tmpFileName);
					std::unique_ptr<image<rgb> > imgConv = LoadImageGenericRGB(tmpFileName);
					double psnr = calculatePSNR(imgIn.get(), imgConv.get());
					imgConv.reset();
					std::filesystem::remove(tmpFileName);
					double bpp = static_cast<double>(8ULL * jpegSize) / static_cast<double>(imgIn->width() * imgIn->height());
					std::cout << std::format("JPEG Quality {} PSNR {} bpp {} total KB {}", quality, psnr, bpp, jpegSize / 1024ULL) << std::endl;
					fileOut << std::format("{}, {}, {}, {}, {}, {}", file, "JPEG", quality, jpegSize, bpp, psnr) << std::endl;
				}
				for (double bppAllocation = 8.0; bppAllocation >= 1.0; bppAllocation -= 1.0) {
					std::unique_ptr<CompressionContextFast> context = createCompressionContextFast(K, BlockSize, bppAllocation);
					std::cout << std::format("MP processing: {} {}%", file, (100L * count) / files.size()) << std::endl;
					size_t outputBytesSize;
					std::unique_ptr<uint8_t[]> encodedBytes = encodeImageFast(
						imgIn.get(),
						K,
						BlockSize,
						context->Y.Quant, context->U.Quant, context->V.Quant,
						context->Y.Dynamic, context->U.Dynamic, context->V.Dynamic,
						outputBytesSize);
					std::unique_ptr<image<rgb> > imgOut = decodeImageFast(encodedBytes.get(), outputBytesSize);
					double psnr = calculatePSNR(imgIn.get(), imgOut.get());
					imgOut.reset();
					double bpp = static_cast<double>(8ULL * outputBytesSize) / static_cast<double>(imgIn->width() * imgIn->height());
					std::cout << std::format("MP Quality {} PSNR {} bpp {} total KB {}", bppAllocation, psnr, bpp, outputBytesSize) << std::endl;
					fileOut << std::format("{}, {}, {}, {}, {}, {}", file, "MP", bppAllocation, outputBytesSize, bpp, psnr) << std::endl;
				}
				imgIn.reset();
			}
			fileOut.flush();
		}
		fileOut.close();
	} else {
		printUsage();
		return -1;
	}
	return 0;
}

