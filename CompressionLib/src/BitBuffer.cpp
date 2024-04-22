#include "../inc/BitBuffer.h"
#include <stdexcept>
#include <cstdlib>
#include <bit>

namespace bitbuffer {

	BitBuffer::BitBuffer() : BitBuffer(32ULL) {
	}

	BitBuffer::BitBuffer(size_t capacity)
		: m_bits(),
		m_writeWordPos(0),
		m_writeBitPos(0),
		m_readWordPos(0),
		m_readBitPos(0)
	{
		m_bits.reserve(capacity);
		m_bits.push_back(0ULL);
		m_bits.push_back(0ULL);
	}

	BitBuffer::BitBuffer(const BitBuffer& other)
		: m_bits(other.m_bits),
		m_writeWordPos(other.m_writeWordPos),
		m_writeBitPos(other.m_writeBitPos),
		m_readWordPos(other.m_readWordPos),
		m_readBitPos(other.m_readBitPos)
	{

	}

	BitBuffer::BitBuffer(BitBuffer&& other) noexcept
		: m_bits(std::move(other.m_bits)),
		m_writeWordPos(other.m_writeWordPos),
		m_writeBitPos(other.m_writeBitPos),
		m_readWordPos(other.m_readWordPos),
		m_readBitPos(other.m_readBitPos)
	{
		other.m_writeWordPos = 0;
		other.m_writeBitPos = 0;
		other.m_readWordPos = 0;
		other.m_readBitPos = 0;
	}

	BitBuffer::~BitBuffer() {
		Clear();
	}

	BitBuffer& BitBuffer::operator=(const BitBuffer& other) {
		if (this != &other)
		{
			m_writeWordPos = other.m_writeWordPos;
			m_writeBitPos = other.m_writeBitPos;
			m_readWordPos = other.m_readWordPos;
			m_readBitPos = other.m_readBitPos;
			m_bits = other.m_bits;
		}
		return *this;
	}

	BitBuffer& BitBuffer::operator=(BitBuffer&& other) noexcept {
		if (this != &other)
		{
			m_writeWordPos = other.m_writeWordPos;
			m_writeBitPos = other.m_writeBitPos;
			m_readWordPos = other.m_readWordPos;
			m_readBitPos = other.m_readBitPos;
			m_bits = std::move(other.m_bits);
			other.m_writeWordPos = 0;
			other.m_writeBitPos = 0;
			other.m_readWordPos = 0;
			other.m_readBitPos = 0;
		}
		return *this;
	}

	void BitBuffer::WriteBits(uint64_t value, int width) {
		if (width < 0 || width > 64) {
			throw new std::range_error("invalid bit width " + width);
		}
		if (width == 0) {
			return;
		}
		if (m_writeWordPos + 1 >= m_bits.size()) {
			m_bits.push_back(0ULL);
			m_bits.push_back(0ULL);
		}
		size_t endpos = m_writeBitPos + width;
		uint64_t mask = (1ULL << width) - 1ULL;
		if (width == 64) {
			mask = 0xFFFFFFFFFFFFFFFFULL;
		}
		uint64_t maskedValue = (value & mask);
		if (endpos <= 64) {
			uint64_t shiftedValue = maskedValue << (64 - endpos);
			uint64_t slot = m_bits[m_writeWordPos];
			m_bits[m_writeWordPos] = slot | shiftedValue;
		} else {
			uint64_t shiftedValue1 = maskedValue >> (endpos - 64);
			uint64_t slot1 = m_bits[m_writeWordPos];
			m_bits[m_writeWordPos] = slot1 | shiftedValue1;
			m_bits[m_writeWordPos + 1] = maskedValue << (128 - endpos);
		}
		m_writeBitPos += width;
		if (m_writeBitPos >= 64) {
			m_writeBitPos -= 64;
			++m_writeWordPos;
		}
	}

	uint64_t BitBuffer::ReadBits(int width) {
		uint64_t retval = PeekBits(width);
		SkipBits(width);
		return retval;
	}

	uint64_t BitBuffer::PeekBits(int width) {
		if (width < 0 || width > 64) {
			throw new std::range_error("invalid bit width " + width);
		}
		width = std::min(width, static_cast<int>((m_writeWordPos - m_readWordPos) * 64ULL + m_writeBitPos - m_readBitPos));
		if (width == 0) {
			return 0;
		}
		size_t endpos = m_readBitPos + width;
		uint64_t mask = (1ULL << width) - 1ULL;
		if (width == 64) {
			mask = 0xFFFFFFFFFFFFFFFFULL;
		}
		uint64_t retval;
		if (endpos <= 64) {
			retval = (m_bits[m_readWordPos] >> (64 - endpos)) & mask;
		} else {
			uint64_t part1 = m_bits[m_readWordPos] << (endpos - 64);
			uint64_t part2 = m_bits[m_readWordPos + 1] >> (128 - endpos);
			retval = (part1 | part2) & mask;
		}
		return retval;
	}

	void BitBuffer::SkipBits(int width) {
		if (width < 0 || width > 64) {
			throw new std::range_error("invalid bit width " + width);
		}
		width = std::min(width, static_cast<int>((m_writeWordPos - m_readWordPos) * 64ULL + m_writeBitPos - m_readBitPos));
		if (width == 0) {
			return;
		}
		m_readBitPos += width;
		if (m_readBitPos >= 64) {
			m_readBitPos -= 64;
			++m_readWordPos;
		}
	}


	void BitBuffer::ResetReadPtr() {
		m_readWordPos = 0;
		m_readBitPos = 0;
	}

	void BitBuffer::Clear() {
		m_bits.clear();
		m_writeWordPos = 0;
		m_writeBitPos = 0;
		m_readWordPos = 0;
		m_readBitPos = 0;
	}

	void BitBuffer::Append(BitBuffer& source) {
		size_t remaining = source.Remaining();
		while (remaining > 0ULL) {
			if (remaining >= 64ULL) {
				WriteLong(source.ReadLong());
			} else {
				WriteBits(source.ReadBits(static_cast<int32_t>(remaining)), static_cast<int32_t>(remaining));
			}
			remaining = source.Remaining();
		}
	}

	std::unique_ptr<uint8_t[]> BitBuffer::Save(size_t& byteLength) const {
		byteLength = 8ULL * m_writeWordPos + ((m_writeBitPos + 7ULL) / 8ULL);
		std::unique_ptr<uint8_t[]> data = std::unique_ptr<uint8_t[]>(new uint8_t[byteLength]);
		for (size_t offset = 0; offset < m_writeWordPos; ++offset) {
			uint64_t value = m_bits[offset];
			data[8ULL * offset] = static_cast<uint8_t>(value >> 56);
			data[8ULL * offset + 1] = static_cast<uint8_t>(value >> 48);
			data[8ULL * offset + 2] = static_cast<uint8_t>(value >> 40);
			data[8ULL * offset + 3] = static_cast<uint8_t>(value >> 32);
			data[8ULL * offset + 4] = static_cast<uint8_t>(value >> 24);
			data[8ULL * offset + 5] = static_cast<uint8_t>(value >> 16);
			data[8ULL * offset + 6] = static_cast<uint8_t>(value >> 8);
			data[8ULL * offset + 7] = static_cast<uint8_t>(value);
		}
		uint64_t lastPart = m_bits[m_writeWordPos];
		for (size_t partOffset = 0; partOffset < ((m_writeBitPos + 7ULL) / 8ULL); ++partOffset) {
			data[8ULL * m_writeWordPos + partOffset] = static_cast<uint8_t>(lastPart >> (8 * (7 - partOffset)));
		}
		return data;
	}

	void BitBuffer::Load(const uint8_t bytes[], size_t startByteOffset, size_t bitLength) {
		m_bits.resize(1 + ((bitLength + 63) / 64));
		m_readWordPos = 0;
		m_readBitPos = 0;
		m_writeWordPos = bitLength / 64ULL;
		m_writeBitPos = bitLength % 64ULL;
		for (size_t offset = 0; offset < m_writeWordPos; ++offset) {
			uint64_t value = (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset]) & 0xFFULL) << 56;
			value |= (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset + 1]) & 0xFFULL) << 48;
			value |= (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset + 2]) & 0xFFULL) << 40;
			value |= (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset + 3]) & 0xFFULL) << 32;
			value |= (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset + 4]) & 0xFFULL) << 24;
			value |= (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset + 5]) & 0xFFULL) << 16;
			value |= (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset + 6]) & 0xFFULL) << 8;
			value |= (static_cast<uint64_t>(bytes[startByteOffset + 8ULL * offset + 7]) & 0xFFULL);
			m_bits[offset] = value;
		}
		uint64_t last = 0ULL;
		for (size_t residual = 0; residual < ((m_writeBitPos + 7ULL) / 8ULL); ++residual) {
			last = last | ((static_cast<uint64_t>(bytes[startByteOffset + (8ULL * m_writeWordPos) + residual]) & 0xFFULL) << (8 * (7 - residual)));
		}
		m_bits[m_writeWordPos] = last;
	}

	void writeGolombCode(uint32_t value,uint32_t M, BitBuffer& buffer) {
		std::div_t parts = std::div(static_cast<int32_t>(value), static_cast<int32_t>(M));
		for (uint32_t i = 0; i < static_cast<uint32_t>(parts.quot); ++i) {
			buffer.WriteBits(1, 1);
		}
		buffer.WriteBits(0, 1);
		uint32_t b = std::bit_width(M);
		uint32_t limit = (1 << (b + 1)) - M;
		if (static_cast<uint32_t>(parts.rem) < limit) {
			buffer.WriteBits(parts.rem, b);
		} else {
			buffer.WriteBits(static_cast<uint64_t>(parts.rem) + limit, b + 1);
		}
	}

	uint32_t readGolombCode(uint32_t M, BitBuffer& buffer) {
		uint32_t q = 0;
		while (buffer.ReadBits(1) != 0) {
			++q;
		}
		uint32_t b = std::bit_width(M);
		uint32_t limit = (1 << (b + 1)) - M;
		uint32_t rem1 = static_cast<uint32_t>(buffer.ReadBits(b));
		uint32_t rem;
		if (rem1 < limit) {
			rem = rem1;
		} else {
			rem = (rem1 << 1) + static_cast<uint32_t>(buffer.ReadBits(1)) - limit;
		}
		return q * M + rem;
	}

	uint32_t golombCodeLength(uint32_t value, uint32_t M) {
		std::div_t parts = std::div(static_cast<int32_t>(value), static_cast<int32_t>(M));
		uint32_t b = std::bit_width(M);
		uint32_t limit = (1 << (b + 1)) - M;
		if (static_cast<uint32_t>(parts.rem) < limit) {
			return b + parts.quot + 1;
		} else {
			return b + parts.quot + 2;
		}
	}

	void writeEliasCode(uint32_t value, BitBuffer& buffer) {
		uint32_t shifted = value + 1;
		uint32_t b = std::bit_width(shifted);
		buffer.WriteBits(0, b - 1);
		buffer.WriteBits(1, 1);
		buffer.WriteBits(shifted, b - 1);
	}

	uint32_t readEliasCode(BitBuffer& buffer) {
		uint32_t n = 0;
		while (buffer.ReadBits(1) == 0) {
			++n;
		}
		return ((1U << n) | static_cast<uint32_t>(buffer.ReadBits(n))) - 1;
	}

	uint32_t eliasCodeLength(uint32_t value) {
		uint32_t b = std::bit_width(value + 1);
		return 2 * b - 1;
	}

	void writeEliasFanoSequenceCode(const std::vector<uint16_t>& nonDecreasingSequence, uint16_t maxSymbol, BitBuffer& buffer) {
		if (nonDecreasingSequence.empty()) {
			return;
		}
		uint32_t m = std::bit_width(maxSymbol);
		uint32_t n = static_cast<uint32_t>(nonDecreasingSequence.size());
		uint32_t nb = std::bit_width(n);
		uint32_t lb = (m >= nb) ? m - nb : 0;
		uint16_t prev = nonDecreasingSequence.front();
		uint32_t prevBucket = 0;
		for (const uint16_t sym : nonDecreasingSequence) {
			if (sym < prev) {
				throw new std::range_error("Sequence must not contain any decreasing values");
			}
			prev = sym;
			uint16_t bucket = sym >> lb;
			while (prevBucket != bucket) {
				buffer.WriteBits(0, 1);
				++prevBucket;
			}
			buffer.WriteBits(1, 1);
		}
		for (const uint16_t sym : nonDecreasingSequence) {
			buffer.WriteBits(sym, lb);
		}
	}

	std::vector<uint16_t> readEliasFanoSequenceCode(size_t sequenceLength, uint16_t maxSymbol, BitBuffer& buffer) {
		std::vector<uint16_t> result(sequenceLength);
		if (sequenceLength == 0) {
			return result;
		}
		uint32_t m = std::bit_width(maxSymbol);
		uint32_t n = static_cast<uint32_t>(sequenceLength);
		uint32_t nb = std::bit_width(n);
		uint32_t lb = (m >= nb) ? m - nb : 0;
		uint16_t bucket = 0;
		for (uint16_t& sym: result) {
			if (buffer.Remaining() == 0ULL) {
				throw new std::range_error("Invalid bitstream");
			}
			while (buffer.ReadBits(1) == 0) {
				++bucket;
			}
			sym = bucket << lb;
		}
		for (uint16_t& sym : result) {
			sym |= static_cast<uint16_t>(buffer.ReadBits(lb));
		}

		return result;
	}

	uint32_t eliasFanoSequenceCodeLength(size_t sequenceLength, uint16_t maxSymbol) {
		if (sequenceLength == 0ULL) {
			return 0;
		}
		uint32_t m = std::bit_width(maxSymbol);
		uint32_t n = static_cast<uint32_t>(sequenceLength);
		uint32_t nb = std::bit_width(n);
		uint32_t lb = (m >= nb) ? m - nb : 0;
		return n * (1 + lb) + (1U << (m - lb)) - 1;
	}
}