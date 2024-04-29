#ifndef _HUFFMAN_H_
#define _HUFFMAN_H_

#include "framework.h"

#include "BitBuffer.h"
#include <memory>
#include <vector>

namespace huffman {

	std::vector<uint16_t> runLengthEncode(const std::vector<uint16_t>& data);
	std::vector<uint16_t> runLengthDecode(const std::vector<uint16_t>& data);

	std::unique_ptr<uint8_t[]> huffmanEncode(const std::vector<uint16_t>& data, size_t& outputLength);
	void huffmanEncode(const std::vector<uint16_t>& data, bitbuffer::BitBuffer& buffer);

	std::vector<uint16_t> huffmanDecode(const uint8_t encoded[], size_t compressedLength);
	void huffmanDecode(bitbuffer::BitBuffer& buffer, std::vector<uint16_t>& decoded);
}

#endif //_HUFFMAN_H_