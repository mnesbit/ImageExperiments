#include "../inc/CompressedImage.h"
#include "../inc/BitBuffer.h"
#include "../inc/Huffman.h"
#include "../../ImageHelper/inc/misc.h"

using namespace bitbuffer;
using namespace img;
using namespace matching;

namespace compressed {

	const uint32_t MAGIC = 0x34317877;

	std::unique_ptr<uint8_t[]> writeCompressed(
		const size_t K,
		const size_t width,
		const size_t height,
		std::vector<uint16_t>& lengths,
		const std::vector<std::vector<uint16_t> >& codes,
		size_t& outputByteSize)
	{
		BitBuffer bitsOut;
		bitsOut.WriteInt(MAGIC);
		bitsOut.WriteInt(static_cast<uint32_t>(width));
		bitsOut.WriteInt(static_cast<uint32_t>(height));
		huffman::huffmanEncode(lengths, bitsOut);
		for (size_t i = 0; i < 2ULL * 3ULL * K; ++i) {
			bitbuffer::BitBuffer tempBits;
			huffman::huffmanEncode(codes[i], tempBits);
			size_t best = tempBits.Remaining();
			int bestFit = -1;
			int m = 1;
			while (m < 2048) {
				size_t estimate = 16ULL;
				for (const auto& sym : codes[i]) {
					estimate += golombCodeLength(static_cast<uint32_t>(sym), m);
				}
				if (estimate < best) {
					best = estimate;
					bestFit = m;
				}
				if ((m & 1) == 1) {
					++m;
				} else {
					m = (m << 1) - 1;
				}
			}
			if (bestFit == -1) {
				bitsOut.WriteBits(0, 1);
				bitsOut.Append(tempBits);
			} else {
				bitsOut.WriteBits(1, 1);
				bitsOut.WriteShort(bestFit);
				for (const auto& sym : codes[i]) {
					writeGolombCode(static_cast<uint32_t>(sym), bestFit, bitsOut);
				}
			}
		}
		return bitsOut.Save(outputByteSize);
	}

	std::unique_ptr<uint8_t[]> encodeImage(const image<rgb>* imgIn, const size_t K, const size_t blockSize,
		const double quantY[], const double quantU[], const double quantV[],
		DynamicDictionaryFunction dynamicY, DynamicDictionaryFunction dynamicU, DynamicDictionaryFunction dynamicV,
		size_t& outputByteSize) {
		const size_t width = imgIn->width();
		const size_t height = imgIn->height();
		std::vector<BasisChoice> choicesY(K);
		std::vector<BasisChoice> choicesU(K);
		std::vector<BasisChoice> choicesV(K);
		math::Vector blockY(blockSize * blockSize);
		math::Vector blockU(blockSize * blockSize);
		math::Vector blockV(blockSize * blockSize);
		std::vector<uint16_t> lengths;
		std::vector<std::vector<uint16_t> > codes(2 * 3 * K);
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
				int countY = CalcMPDynamic(static_cast<int>(K), quantY, choicesY, blockY, dynamicY);
				lengths.push_back(countY);
				for (int i = 0; i < countY; ++i) {
					codes[2ULL * i].push_back(choicesY[i].deltaId);
					codes[2ULL * i + 1ULL].push_back(choicesY[i].intCoeff);
				}
				int countU = CalcMPDynamic(static_cast<int>(K), quantU, choicesU, blockU, dynamicU);
				lengths.push_back(countU);
				for (int i = 0; i < countU; ++i) {
					codes[2ULL * K + 2ULL * i].push_back(choicesU[i].deltaId);
					codes[2ULL * K + 2ULL * i + 1ULL].push_back(choicesU[i].intCoeff);
				}
				int countV = CalcMPDynamic(static_cast<int>(K), quantV, choicesV, blockV, dynamicV);
				lengths.push_back(countV);
				for (int i = 0; i < countV; ++i) {
					codes[4ULL * K + 2ULL * i].push_back(choicesV[i].deltaId);
					codes[4ULL * K + 2ULL * i + 1ULL].push_back(choicesV[i].intCoeff);
				}
			}
		}
		return writeCompressed(K, width, height, lengths, codes, outputByteSize);
	}

	void readCompressed(
		const uint8_t bytes[],
		size_t byteSize,
		const size_t K,
		size_t& width,
		size_t& height,
		std::vector<uint16_t>& lengths,
		std::vector<std::vector<uint16_t> >& codes
	) {
		BitBuffer bitsIn;
		bitsIn.Load(bytes, 0ULL, 8ULL * byteSize);
		if (bitsIn.ReadInt() != MAGIC) {
			throw new std::range_error("Invalid input data");
		}
		width = static_cast<size_t>(bitsIn.ReadInt());
		height = static_cast<size_t>(bitsIn.ReadInt());
		huffman::huffmanDecode(bitsIn, lengths);
		for (size_t i = 0; i < 2ULL * 3ULL * K; ++i) {
			if (bitsIn.ReadBits(1) == 0) {
				huffman::huffmanDecode(bitsIn, codes[i]);
			} else {
				size_t layer = (i / 2ULL) / K;
				size_t depth = (i / 2ULL) % K;
				int m = bitsIn.ReadShort();
				for (size_t offset = 0; offset < (lengths.size() / 3ULL); ++offset) {
					if (lengths[(3ULL * offset) + layer] > static_cast<uint16_t>(depth)) {
						codes[i].push_back(static_cast<uint16_t>(readGolombCode(m, bitsIn)));
					}
				}
			}
		}
	}

	image<rgb>* decodeImage(const size_t K, const size_t blockSize, const uint8_t bytes[], size_t byteSize,
		const double quantY[], const double quantU[], const double quantV[],
		DynamicDictionaryFunction dynamicY, DynamicDictionaryFunction dynamicU, DynamicDictionaryFunction dynamicV) {
		std::vector<uint16_t> lengths;
		std::vector<std::vector<uint16_t> > codes(2 * 3 * K);
		size_t width;
		size_t height;
		readCompressed(bytes, byteSize, K, width, height, lengths, codes);
		image<rgb>* imgOut = new image<rgb>(width, height, false);
		std::vector<size_t> offsets(3 * K);
		size_t lengthOffset = 0;
		std::vector<BasisChoice> choicesY(K);
		std::vector<BasisChoice> choicesU(K);
		std::vector<BasisChoice> choicesV(K);
		for (size_t x = 0; x < width; x += blockSize) {
			std::cout << std::format("decode {} %", (100 * x) / width) << std::endl;
			for (size_t y = 0; y < height; y += blockSize) {
				int countY = lengths[lengthOffset++];
				for (int i = 0; i < countY; ++i) {
					choicesY[i].deltaId = codes[2ULL * i][offsets[i]];
					choicesY[i].intCoeff = codes[2ULL * i + 1ULL][offsets[i]++];
				}
				math::Vector decodedY = FromCoeffsDynamic(countY, quantY, choicesY, dynamicY);
				int countU = lengths[lengthOffset++];
				for (int i = 0; i < countU; ++i) {
					choicesU[i].deltaId = codes[2ULL * K + 2ULL * i][offsets[K + i]];
					choicesU[i].intCoeff = codes[2ULL * K + 2ULL * i + 1ULL][offsets[K + i]++];
				}
				math::Vector decodedU = FromCoeffsDynamic(countU, quantU, choicesU, dynamicU);
				int countV = lengths[lengthOffset++];
				for (int i = 0; i < countV; ++i) {
					choicesV[i].deltaId = codes[4ULL * K + 2ULL * i][offsets[2ULL * K + i]];
					choicesV[i].intCoeff = codes[4ULL * K + 2ULL * i + 1ULL][offsets[2ULL * K + i]++];
				}
				math::Vector decodedV = FromCoeffsDynamic(countV, quantV, choicesV, dynamicV);
				for (size_t dx = 0; dx < blockSize; ++dx) {
					size_t u = x + dx;
					for (size_t dy = 0; dy < blockSize; ++dy) {
						size_t v = y + dy;
						if ((u < width) && (v < height)) {
							rgb pt = RGBFromYUV(yuv{
								.y = decodedY[dx + (blockSize * dy)],
								.u = decodedU[dx + (blockSize * dy)],
								.v = decodedV[dx + (blockSize * dy)]
								});
							imRef(imgOut, u, v) = pt;
						}
					}
				}
			}
		}
		return imgOut;
	}

	image<rgb>* decodeImageWithError(const size_t K, const size_t blockSize, const uint8_t bytes[], size_t byteSize,
		const double quantY[], const double quantU[], const double quantV[],
		DynamicDictionaryFunction dynamicY, DynamicDictionaryFunction dynamicU, DynamicDictionaryFunction dynamicV,
		const image<rgb>* original, double& squaredError) {
		image<rgb>* imgOut = decodeImage(K, blockSize, bytes, byteSize, quantY, quantU, quantV, dynamicY, dynamicU, dynamicV);
		for (size_t x = 0; x < imgOut->width(); ++x) {
			for (size_t y = 0; y < imgOut->height(); ++y) {
				rgb pt = imRef(imgOut, x, y);
				rgb ptOrig = imRef(original, x, y);
				double errR = static_cast<double>(ptOrig.r) - static_cast<double>(pt.r);
				double errG = static_cast<double>(ptOrig.g) - static_cast<double>(pt.g);
				double errB = static_cast<double>(ptOrig.b) - static_cast<double>(pt.b);
				squaredError += square(errR) + square(errG) + square(errB);
			}
		}
		return imgOut;
	}

}