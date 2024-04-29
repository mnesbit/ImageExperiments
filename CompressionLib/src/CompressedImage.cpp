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

	const double Y_DECAY = 0.903273474567349;
	static const double s_varY[32] = {
		1449064.59233812,
		30820.1579762722,
		4844.63260676556,
		2053.92773816285,
		1172.9896493733,
		753.237119689743,
		518.819397824187,
		375.240182302887,
		281.147739274468,
		216.685057541455,
		170.900839234616,
		137.390795170594,
		112.155972024385,
		92.9252736923566,
		77.8760164175741,
		65.9459337614281,
		56.3353591922793,
		48.4368336883129,
		41.925891589094,
		36.4751772230096,
		31.8461765285566,
		27.8910236991453,
		24.5017939480865,
		21.565734506665,
		19.0155226299456,
		16.7729497343868,
		14.8077163462738,
		13.0709717132486,
		11.5415673988584,
		10.181296913059,
		8.98377672266684,
		7.91893798608267
	};

	const double U_DECAY = 0.896203505769308;
	static const double s_varU[32] = {
		52077.849444871,
		817.230961962173,
		108.460079765957,
		49.6660006852308,
		28.7579729013414,
		19.079159877342,
		13.7299518406063,
		10.4400924202154,
		8.2619144137916,
		6.7267078724334,
		5.57422602145882,
		4.7025800480553,
		4.01946465126433,
		3.46931478622891,
		3.01828938492078,
		2.63983764986421,
		2.32368370313153,
		2.05390567923774,
		1.8238006064889,
		1.62570297252752,
		1.45113342100199,
		1.29836265266357,
		1.16390680518714,
		1.04498557529108,
		0.942864343639315,
		0.850429307744034,
		0.767969071651754,
		0.697079095720249,
		0.632642685143598,
		0.577207756783642,
		0.527603303249553,
		0.484287449680384
	};

	const double V_DECAY = 0.894385597563333;
	static const double s_varV[32] = {
		60533.3476624072,
		616.171423884389,
		70.4885137403371,
		31.192929256634,
		16.7701486843641,
		10.4440325431894,
		7.14761351395262,
		5.23282665460927,
		4.01950511465798,
		3.21953682878113,
		2.65809302662582,
		2.25228838110326,
		1.95339412948984,
		1.72246574320484,
		1.54043659078564,
		1.3955148561969,
		1.27896502556882,
		1.17630144588298,
		1.09338029569913,
		1.01735307088158,
		0.951127331040105,
		0.89334028890154,
		0.839607375373023,
		0.792835347361864,
		0.744840989369572,
		0.704372940300665,
		0.667291886802588,
		0.630696739701261,
		0.596960624803531,
		0.567088882662493,
		0.538446158863587,
		0.51091395549458
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