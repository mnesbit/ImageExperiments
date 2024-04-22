#ifndef _COMPRESSED_IMAGE_H_
#define _COMPRESSED_IMAGE_H_

#include <memory>
#include "MatchingPursuit.h"
#include "../../ImageHelper/inc/image.h"

namespace compressed {
	std::unique_ptr<uint8_t[]> encodeImage(
		const img::image<img::rgb>* imgIn,
		const size_t K,
		const size_t blockSize,
		const double quantY[], const double quantU[], const double quantV[],
		matching::DynamicDictionaryFunction dynamicY, matching::DynamicDictionaryFunction dynamicU, matching::DynamicDictionaryFunction dynamicV,
		size_t& outputByteSize);

	img::image<img::rgb>* decodeImage(
		const size_t K,
		const size_t blockSize,
		const uint8_t bytes[],
		size_t byteSize,
		const double quantY[], const double quantU[], const double quantV[],
		matching::DynamicDictionaryFunction dynamicY, matching::DynamicDictionaryFunction dynamicU, matching::DynamicDictionaryFunction dynamicV);

	img::image<img::rgb>* decodeImageWithError(
		const size_t K,
		const size_t blockSize,
		const uint8_t bytes[],
		size_t byteSize,
		const double quantY[], const double quantU[], const double quantV[],
		matching::DynamicDictionaryFunction dynamicY, matching::DynamicDictionaryFunction dynamicU, matching::DynamicDictionaryFunction dynamicV,
		const img::image<img::rgb>* original, double& squaredError);
}

#endif //_COMPRESSED_IMAGE_H_