#ifndef _BIT_BUFFER_H_
#define _BIT_BUFFER_H_

#include <vector>
#include <iostream>
#include <bitset>
#include <format>

namespace bitbuffer {

	inline int32_t zigzagDecode(uint32_t x) {
		return static_cast<int32_t>(((x >> 1) ^ -(static_cast<int64_t>(x) & 1)));
	}

	inline uint32_t zigzagEncode(int32_t x) {
		return (static_cast<uint32_t>(x) << 1) ^ (x >> 31);
	}

	class BitBuffer {
	public:
		BitBuffer();
		BitBuffer(size_t capacity);
		BitBuffer(const BitBuffer& other);
		BitBuffer(BitBuffer&& other) noexcept;
		~BitBuffer();

		BitBuffer& operator=(const BitBuffer& other);
		BitBuffer& operator=(BitBuffer&& other) noexcept;

		size_t Size() const {
			return 64ULL * m_writeWordPos + m_writeBitPos;
		}

		size_t ByteSize() const {
			return ((Size() + 7ULL) / 8ULL);
		}

		void WriteBits(uint64_t value, int width);
		uint64_t ReadBits(int width);
		uint64_t PeekBits(int width);
		void SkipBits(int width);

		void WriteByte(uint8_t value) {
			WriteBits(static_cast<uint64_t>(value), 8);
		}

		uint8_t ReadByte() {
			return static_cast<uint8_t>(ReadBits(8));
		}

		void WriteShort(uint16_t value) {
			WriteBits(static_cast<uint64_t>(value), 16);
		}

		uint16_t ReadShort() {
			return static_cast<uint16_t>(ReadBits(16));
		}

		void WriteInt(uint32_t value) {
			WriteBits(static_cast<uint64_t>(value), 32);
		}

		uint32_t ReadInt() {
			return static_cast<uint32_t>(ReadBits(32));
		}

		void WriteLong(uint64_t value) {
			WriteBits(value, 64);
		}

		uint64_t ReadLong() {
			return ReadBits(64);
		}

		std::unique_ptr<uint8_t[]> Save(size_t& byteLength) const;
		void Load(const uint8_t bytes[], size_t startByteOffset, size_t bitLength);

		void ResetReadPtr();
		void Clear();


		friend std::ostream& operator<<(std::ostream& os, BitBuffer const& bits) {
			for (size_t offset = 0; offset <= bits.m_writeWordPos; ++offset) {
				os << std::format("{:016X}", bits.m_bits[offset]) << std::endl;
			}
			//for (size_t offset = 0; offset < bits.m_writeWordPos; ++offset) {
			//	std::bitset<64> slot(bits.m_bits[offset]);
			//	for (size_t bitpos = 0; bitpos < 64; ++bitpos) {
			//		os << slot[bitpos];
			//	}
			//	os << std::endl;
			//}
			//std::bitset<64> tail(bits.m_bits[bits.m_writeWordPos]);
			//for (size_t bitpos = 0; bitpos < bits.m_writeBitPos; ++bitpos) {
			//	os << tail[bitpos];
			//}
			//os << std::endl;
			return os;
		}
	private:
		std::vector<uint64_t> m_bits;
		size_t m_writeWordPos;
		size_t m_writeBitPos;
		size_t m_readWordPos;
		size_t m_readBitPos;
	};

	void writeGolombCode(uint32_t value, uint32_t M, BitBuffer& buffer);
	uint32_t readGolombCode(uint32_t M, BitBuffer& buffer);
	uint8_t golombCodeLength(uint32_t value, uint32_t M);

	void writeEliasCode(uint32_t value, BitBuffer& buffer);
	uint32_t readEliasCode(BitBuffer& buffer);
	uint8_t eliasCodeLength(uint32_t value);

}

#endif// _BIT_BUFFER_H_