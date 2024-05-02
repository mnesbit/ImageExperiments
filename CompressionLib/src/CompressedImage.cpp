#include "../inc/CompressedImage.h"
#include "../inc/BitBuffer.h"
#include "../inc/Huffman.h"
#include "../inc/BasisSet.h"
#include "../../ImageHelper/inc/misc.h"

using namespace bitbuffer;
using namespace img;
using namespace matching;
using namespace basis;

namespace compressed {

	const uint32_t MAGIC = 0x4D4E3234;

	const double Y_DECAY = 0.902045039488061;
	static const double s_varY[32] = {
		1449455.61399403,
		30867.8722232759,
		4879.76236869648,
		2065.81004100418,
		1177.78544096912,
		754.545827240537,
		519.229145237154,
		375.509017928094,
		281.02055698585,
		216.291896608802,
		170.433377219481,
		137.390795170594,
		111.760859784514,
		92.465892223227,
		77.4762763021059,
		65.5657474800836,
		55.9694631412915,
		48.1296222409122,
		41.6533463654379,
		36.2103188098409,
		31.620636785048,
		27.6994619197188,
		24.3350442530727,
		21.416843198437,
		18.879120562683,
		16.6557974119762,
		14.7004164720848,
		12.9791457307132,
		11.4562911334688,
		10.1096285629269,
		8.9223405657796,
		7.86882234976579
	};

	const double U_DECAY = 0.896332644824969;
	static const double s_varU[32] = {
		51995.6231219068,
		814.727839313831,
		108.677634702502,
		49.7911952400269,
		28.7429354781437,
		19.0770865041308,
		13.7134473152652,
		10.4220864760027,
		8.24661056024952,
		6.70462984374069,
		5.55959248319828,
		4.6929026981696,
		4.00644915447918,
		3.4569689589025,
		3.00626981236079,
		2.63088816527246,
		2.31483558804887,
		2.04636560444097,
		1.81780225821739,
		1.61760270635648,
		1.44182998235134,
		1.29183839884355,
		1.15688638919911,
		1.0404065817416,
		0.938042133945239,
		0.844584472128933,
		0.764573763603021,
		0.691885170456949,
		0.629607793503852,
		0.573021503073893,
		0.52497651127712,
		0.483288166342945
	};

	const double V_DECAY = 0.897340618787505;
	static const double s_varV[32] = {
		60578.6241767756,
		617.61939120778,
		70.9553277465942,
		31.4166652349442,
		16.8206825627114,
		10.4578171000126,
		7.14592512982323,
		5.22663706003167,
		4.01504652955091,
		3.20518029980816,
		2.64504709802794,
		2.2477896896281,
		1.94630992425302,
		1.72002308826788,
		1.53680040297081,
		1.39452712284538,
		1.27393363160348,
		1.17553215497528,
		1.08841824408173,
		1.01505869329656,
		0.950595903582374,
		0.893571526043102,
		0.841836849604898,
		0.792090537556817,
		0.74756288452653,
		0.704561007977196,
		0.665399335594415,
		0.631274499723472,
		0.597331898015673,
		0.568852846586831,
		0.539371538069597,
		0.513182721162335
	};

	void createQuantizationTables(
		const size_t K,
		const double bppAllocation,
		math::Vector& quantY, math::Vector& quantU, math::Vector& quantV
	) {
		double allocated = 0.0;
		std::vector<double> alloc(3ULL * K);
		std::vector<double> variances(3ULL * K);
		for (size_t i = 0; i < K; ++i) {
			variances[i] = s_varY[i];
			variances[i + K] = s_varU[i];
			variances[i + 2ULL * K] = s_varV[i];
		}
		while ((allocated / (8.0 * 8.0)) < bppAllocation) {
			auto maxElement = std::max_element(variances.cbegin(), variances.cend());
			ptrdiff_t index = maxElement - variances.cbegin();
			alloc[index] += 1.0;
			variances[index] /= 2.0;
			allocated += 1.0;
		}
		quantY = math::Vector(K);
		quantU = math::Vector(K);
		quantV = math::Vector(K);
		for (size_t i = 0; i < K; ++i) {
			double min = 1.0;
			if (i == 0) {
				min = 8.0;
			}
			quantY[i] = std::max(std::ceil(255.0 * 8.0 * std::pow(Y_DECAY, static_cast<double>(i)) * std::pow(0.5, alloc[i])), min);
			quantU[i] = std::max(std::ceil(255.0 * 8.0 * std::pow(U_DECAY, static_cast<double>(i)) * std::pow(0.5, alloc[i + K])), min);
			quantV[i] = std::max(std::ceil(255.0 * 8.0 * std::pow(V_DECAY, static_cast<double>(i)) * std::pow(0.5, alloc[i + 2ULL * K])), min);
		}
	}

	static math::Matrix dynamicBasis(
		const math::Matrix& baseDict,
		const std::vector<math::Matrix>& detailBasis,
		int prevCoeffs,
		const std::vector<BasisChoice>& prevChoices
	)  {
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
			math::Matrix dictionary(atomCount, baseDict.Columns());
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
					memcpy(pDict + (offset * dictionary.Columns()), detail.Data(), sizeof(double) * detail.Rows() * detail.Columns());
					offset += static_cast<int>(detail.Rows());
				}
			}
			return dictionary;
		};

	std::unique_ptr<CompressionContext> createCompressionContext(size_t K, size_t blockSize, double bppAllocation) {
		std::unique_ptr<CompressionContext> context = std::unique_ptr<CompressionContext>(new CompressionContext);
		context->K = K;
		context->BlockSize = blockSize;
		createQuantizationTables(K, bppAllocation, context->Y.Quant, context->U.Quant, context->V.Quant);
		std::vector<Line> basisSet = distinctLineShapes(blockSize);
		context->BaseDict = createSegmentDictionary(blockSize, basisSet);
		size_t fullRowCount = context->BaseDict.Rows();
		for (const Line& basis : basisSet) {
			context->Y.DetailBasis.push_back(createIntraSegmentDictionary(blockSize, basis, covariancePredictionModelY));
			context->U.DetailBasis.push_back(createIntraSegmentDictionary(blockSize, basis, covariancePredictionModelU));
			context->V.DetailBasis.push_back(createIntraSegmentDictionary(blockSize, basis, covariancePredictionModelV));
		}
		context->Y.Dynamic = [&context](int prevCoeffs, const std::vector<BasisChoice>& prevChoices)  -> math::Matrix {
			return dynamicBasis(context->BaseDict, context->Y.DetailBasis, prevCoeffs, prevChoices);
			};
		context->U.Dynamic = [&context](int prevCoeffs, const std::vector<BasisChoice>& prevChoices)  -> math::Matrix {
			return dynamicBasis(context->BaseDict, context->U.DetailBasis, prevCoeffs, prevChoices);
			};
		context->V.Dynamic = [&context](int prevCoeffs, const std::vector<BasisChoice>& prevChoices)  -> math::Matrix {
			return dynamicBasis(context->BaseDict, context->V.DetailBasis, prevCoeffs, prevChoices);
			};
		return context;
	}

	double calculatePSNR(const img::image<img::rgb>* original, const img::image<img::rgb>* decoded) {
		double squaredError = 0.0;
		for (size_t x = 0ULL; x < original->width(); ++x) {
			for (size_t y = 0ULL; y < original->height(); ++y) {
				rgb pt = imRef(decoded, x, y);
				rgb ptOrig = imRef(original, x, y);
				double r = static_cast<double>(pt.r) - static_cast<double>(ptOrig.r);
				double g = static_cast<double>(pt.g) - static_cast<double>(ptOrig.g);
				double b = static_cast<double>(pt.b) - static_cast<double>(ptOrig.b);
				squaredError += r * r + g * g + b * b;
			}
		}
		double	mse = squaredError / (original->width() * original->height());
		return 20.0 * log10(3.0 * 255.0) - 10.0 * log10(mse);
	}

	void writeHuffmanOrGolomb(const std::vector<uint16_t>& data, BitBuffer& buffer) {
		bitbuffer::BitBuffer tempBits;
		huffman::huffmanEncode(data, tempBits);
		size_t best = tempBits.Remaining();
		int bestFit = -1;
		int m = 1;
		while (m < 2048) {
			size_t estimate = 16ULL;
			for (const auto& sym : data) {
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
			buffer.WriteBits(0, 1);
			buffer.Append(tempBits);
		} else {
			buffer.WriteBits(1, 1);
			buffer.WriteShort(bestFit);
			for (const auto& sym : data) {
				writeGolombCode(static_cast<uint32_t>(sym), bestFit, buffer);
			}
		}
	}

	void readHuffmanOrGolomb(BitBuffer& buffer, size_t length, std::vector<uint16_t>& data) {
		if (buffer.ReadBits(1) == 0) {
			huffman::huffmanDecode(buffer, data);
		} else {
			int m = buffer.ReadShort();
			for (size_t offset = 0; offset < length; ++offset) {
				data.push_back(static_cast<uint16_t>(readGolombCode(m, buffer)));
			}
		}
	}

	std::unique_ptr<uint8_t[]> writeCompressed(
		const size_t K,
		const size_t blockSize,
		const size_t width,
		const size_t height,
		const double quantY[], const double quantU[], const double quantV[],
		std::vector<uint16_t>& lengths,
		const std::vector<std::vector<uint16_t> >& codes,
		size_t& outputByteSize)
	{
		BitBuffer bitsOut;
		bitsOut.WriteInt(MAGIC);
		bitsOut.WriteInt(static_cast<uint32_t>(width));
		bitsOut.WriteInt(static_cast<uint32_t>(height));
		bitsOut.WriteByte(static_cast<uint8_t>(K));
		bitsOut.WriteByte(static_cast<uint8_t>(blockSize));
		for (size_t i = 0; i < K; ++i) {
			bitsOut.WriteShort(static_cast<uint16_t>(quantY[i]));
		}
		for (size_t i = 0; i < K; ++i) {
			bitsOut.WriteShort(static_cast<uint16_t>(quantU[i]));
		}
		for (size_t i = 0; i < K; ++i) {
			bitsOut.WriteShort(static_cast<uint16_t>(quantV[i]));
		}

		writeHuffmanOrGolomb(lengths, bitsOut);
		for (size_t i = 0; i < 2ULL * 3ULL * K; ++i) {
			std::vector<uint16_t> compCodes = huffman::runLengthEncode(codes[i]);
			if (compCodes.size() + 4 < codes[i].size()) {
				bitsOut.WriteBits(1, 1);
				bitsOut.WriteInt(static_cast<uint32_t>(compCodes.size()));
				writeHuffmanOrGolomb(compCodes, bitsOut);
			} else {
				bitsOut.WriteBits(0, 1);
				writeHuffmanOrGolomb(codes[i], bitsOut);
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
		return writeCompressed(K, blockSize, width, height, quantY, quantU, quantV, lengths, codes, outputByteSize);
	}

	std::unique_ptr<CompressionContext> readCompressed(
		const uint8_t bytes[],
		size_t byteSize,
		size_t& K,
		size_t& blockSize,
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
		K = static_cast<size_t>(bitsIn.ReadByte());
		blockSize = static_cast<size_t>(bitsIn.ReadByte());
		std::unique_ptr<CompressionContext> context = createCompressionContext(K, blockSize, 0.0);
		for (size_t i = 0; i < K; ++i) {
			context->Y.Quant[i] = static_cast<double>(bitsIn.ReadShort());
		}
		for (size_t i = 0; i < K; ++i) {
			context->U.Quant[i] = static_cast<double>(bitsIn.ReadShort());
		}
		for (size_t i = 0; i < K; ++i) {
			context->V.Quant[i] = static_cast<double>(bitsIn.ReadShort());
		}
		codes.clear();
		codes.resize(2ULL * 3ULL * K);
		size_t patchLength = 3ULL * ((width + 7ULL) / 8ULL) * ((height + 7ULL) / 8ULL);
		lengths.reserve(patchLength);
		readHuffmanOrGolomb(bitsIn, patchLength, lengths);
		for (size_t i = 0; i < 2ULL * 3ULL * K; ++i) {
			if (bitsIn.ReadBits(1) == 1ULL) {
				size_t codeLength = static_cast<size_t>(bitsIn.ReadInt());
				std::vector<uint16_t> compCodes;
				compCodes.reserve(codeLength);
				readHuffmanOrGolomb(bitsIn, codeLength, compCodes);
				std::vector<uint16_t> decomp = huffman::runLengthDecode(compCodes);
				codes[i].swap(decomp);
			} else {
				size_t layer = (i / 2ULL) / K;
				size_t depth = (i / 2ULL) % K;
				size_t layerLength = 0ULL;
				for (size_t offset = 0; offset < (lengths.size() / 3ULL); ++offset) {
					if (lengths[(3ULL * offset) + layer] > static_cast<uint16_t>(depth)) {
						++layerLength;
					}
				}
				codes[i].reserve(layerLength);
				readHuffmanOrGolomb(bitsIn, layerLength, codes[i]);
			}
		}
		return context;
	}

	image<rgb>* decodeImage(const uint8_t bytes[], size_t byteSize) {
		size_t width;
		size_t height;
		size_t K;
		size_t blockSize;
		std::vector<uint16_t> lengths;
		std::vector<std::vector<uint16_t> > codes;
		std::unique_ptr<CompressionContext> context = readCompressed(bytes, byteSize, K, blockSize, width, height, lengths, codes);
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
				math::Vector decodedY = FromCoeffsDynamic(countY, context->Y.Quant.Data(), choicesY, context->Y.Dynamic);
				int countU = lengths[lengthOffset++];
				for (int i = 0; i < countU; ++i) {
					choicesU[i].deltaId = codes[2ULL * K + 2ULL * i][offsets[K + i]];
					choicesU[i].intCoeff = codes[2ULL * K + 2ULL * i + 1ULL][offsets[K + i]++];
				}
				math::Vector decodedU = FromCoeffsDynamic(countU, context->U.Quant.Data(), choicesU, context->U.Dynamic);
				int countV = lengths[lengthOffset++];
				for (int i = 0; i < countV; ++i) {
					choicesV[i].deltaId = codes[4ULL * K + 2ULL * i][offsets[2ULL * K + i]];
					choicesV[i].intCoeff = codes[4ULL * K + 2ULL * i + 1ULL][offsets[2ULL * K + i]++];
				}
				math::Vector decodedV = FromCoeffsDynamic(countV, context->V.Quant.Data(), choicesV, context->V.Dynamic);
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

}