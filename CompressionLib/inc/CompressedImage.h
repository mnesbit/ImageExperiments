#ifndef _COMPRESSED_IMAGE_H_
#define _COMPRESSED_IMAGE_H_

#include <memory>
#include "MatchingPursuit.h"
#include "../../ImageHelper/inc/image.h"

namespace compressed {
	// created using bit allocation strategy outlined in:
	// Thakur, V. S., Thakur, K., Gupta, S., & Rao, K. R. (2020). Image?independent optimal non?negative integer bit allocation technique for the DCT?based image transform coders. IET Image Processing, 14(1), 11-24.
	// Coefficients then derived from max possible range (i.e. 255 * 8 = 2040 for 8 x 8 patches on 8-bit YUV data)
	// adjusted for the observed rate of range contraction e.g. 0.9^coefficient index
	// divided by 2^bits allocated from above procedure.
	// The DC component quantization is capped at 8.0 and other coefficients are capped at 1.0 divisor
	void createQuantizationTables(
		const size_t K,
		const size_t blockSize,
		const double bppAllocation,
		math::Vector& quantY, math::Vector& quantU, math::Vector& quantV
		);

	struct ChannelContext {
		math::Vector Quant;
		std::vector<math::Matrix> DetailBasis;
		matching::DynamicDictionaryFunction Dynamic;
	};

	struct CompressionContext {
		size_t K{ 32 };
		size_t BlockSize{ 8 };
		math::Matrix BaseDict;
		ChannelContext Y;
		ChannelContext U;
		ChannelContext V;
	};

	struct ChannelContextFast {
	public:
		Eigen::VectorXf Quant;
		std::vector<Eigen::MatrixXf> DetailBasis;
		matching::DynamicDictionaryFunctionFast Dynamic;
	};

	struct CompressionContextFast {
	public:
		size_t K{ 32 };
		size_t BlockSize{ 8 };
		Eigen::MatrixXf BaseDict;
		ChannelContextFast Y;
		ChannelContextFast U;
		ChannelContextFast V;
	};

	std::unique_ptr<CompressionContext> createCompressionContext(size_t K, size_t blockSize, double bppAllocation);
	std::unique_ptr<CompressionContextFast> createCompressionContextFast(size_t K, size_t blockSize, double bppAllocation);

	double calculatePSNR(const img::image<img::rgb>* original, const img::image<img::rgb>* decoded);

	std::unique_ptr<uint8_t[]> encodeImage(
		const img::image<img::rgb>* imgIn,
		const size_t K,
		const size_t blockSize,
		const double quantY[], const double quantU[], const double quantV[],
		const matching::DynamicDictionaryFunction& dynamicY, const matching::DynamicDictionaryFunction& dynamicU, const matching::DynamicDictionaryFunction& dynamicV,
		size_t& outputByteSize);

	std::unique_ptr<uint8_t[]> encodeImageFast(
		const img::image<img::rgb>* imgIn,
		const size_t K,
		const size_t blockSize,
		const Eigen::VectorXf& quantY, const Eigen::VectorXf& quantU, const Eigen::VectorXf& quantV,
		const matching::DynamicDictionaryFunctionFast& dynamicY, const matching::DynamicDictionaryFunctionFast& dynamicU, const matching::DynamicDictionaryFunctionFast& dynamicV,
		size_t& outputByteSize);

	img::image<img::rgb>* decodeImage(const uint8_t bytes[], size_t byteSize);
	img::image<img::rgb>* decodeImageFast(const uint8_t bytes[], size_t byteSize);
}

#endif //_COMPRESSED_IMAGE_H_