#include "pch.h"

#include "../CompressionLib/inc/BitBuffer.h"
#include <random>

TEST(BitBufferTests, EmptyTests) {
	bitbuffer::BitBuffer empty;
	ASSERT_EQ(empty.Size(), 0ULL) << "Should be no bits";
	ASSERT_EQ(empty.ByteSize(), 0ULL) << "Should be no bytes";
	ASSERT_EQ(empty.ReadPos(), 0ULL) << "Should be at 0 position";
	ASSERT_EQ(empty.Remaining(), 0ULL) << "Should be no bits";
	ASSERT_EQ(empty.ReadByte(), static_cast<uint8_t>(0)) << "Reads beyond end should return 0";
	ASSERT_EQ(empty.ReadShort(), static_cast<uint16_t>(0)) << "Reads beyond end should return 0";
	ASSERT_EQ(empty.ReadInt(), static_cast<uint32_t>(0)) << "Reads beyond end should return 0";
	ASSERT_EQ(empty.ReadLong(), static_cast<uint64_t>(0)) << "Reads beyond end should return 0";
	for (int i = 0; i < 64; ++i) {
		ASSERT_EQ(empty.ReadBits(i), 0ULL) << "Reads beyond end should return 0";
	}
	ASSERT_THROW(empty.ReadBits(-1), std::range_error*);
	ASSERT_THROW(empty.ReadBits(65), std::range_error*);
	ASSERT_EQ(empty.Size(), 0ULL) << "Should be no bits";
	ASSERT_EQ(empty.ReadPos(), 0ULL) << "Should be at 0 position";
	ASSERT_EQ(empty.Remaining(), 0ULL) << "Should be no bits";

	size_t size;
	std::unique_ptr<uint8_t[]> saved = empty.Save(size);
	ASSERT_EQ(size, 0ULL) << "Should be no bytes";
	bitbuffer::BitBuffer empty2;
	empty2.WriteBits(1, 1); //make non-empty
	empty2.Load(saved.get(), 0ULL, 0ULL); //should not throw and should reset everything
	ASSERT_EQ(empty2.ReadPos(), 0ULL) << "Should be at 0 position";
	ASSERT_EQ(empty2.Remaining(), 0ULL) << "Should be no bits";
	ASSERT_EQ(empty2.Size(), 0ULL) << "Should be no bits";
	ASSERT_EQ(empty2.ByteSize(), 0ULL) << "Should be no bytes";
}

TEST(BitBufferTests, SimpleReadWriteTests) {
	bitbuffer::BitBuffer buffer;
	size_t length = 0ULL;
	buffer.WriteBits(1, 0);
	ASSERT_EQ(buffer.Size(), 0ULL) << "Should be no bits";
	ASSERT_EQ(buffer.ByteSize(), 0ULL) << "Should be no bytes";
	ASSERT_EQ(buffer.ReadPos(), 0ULL) << "Should be at 0 position";
	ASSERT_EQ(buffer.Remaining(), 0ULL) << "Should be no bits";
	ASSERT_THROW(buffer.ReadBits(-1), std::range_error*);
	ASSERT_THROW(buffer.ReadBits(65), std::range_error*);

	for (int width = 1; width <= 64; ++width) {
		uint64_t value = 0ULL;
		while (value < (1ULL << width)) {
			length += width;
			uint64_t writeValue = value;
			buffer.WriteBits(value, width);
			ASSERT_EQ(buffer.Size(), length) << "write should affect size";
			ASSERT_EQ(buffer.Remaining(), length) << "remaining should increase";
			ASSERT_EQ(buffer.ReadPos(), 0ULL) << "readpos should stay at 0";
			if ((value & 1ULL) == 0) {
				value |= 1ULL;
			} else {
				value = (value << 1) ^ value;
			}
		}
	}
	size_t byteLen;
	std::unique_ptr<uint8_t[]> saved = buffer.Save(byteLen);
	bitbuffer::BitBuffer buffer2;
	buffer2.Load(saved.get(), 0ULL, buffer.Size());
	size_t readLength = 0ULL;
	size_t maxLength = buffer.Size();
	for (int width = 1; width <= 64; ++width) {
		uint64_t value = 0ULL;
		while (value < (1ULL << width)) {
			uint64_t readValue = buffer.ReadBits(width);
			ASSERT_EQ(readValue, value) << "Readback should be of expected values";
			ASSERT_EQ(buffer.Size(), maxLength) << "Size should not change";
			readLength += width;
			ASSERT_EQ(buffer.ReadPos(), readLength) << "Read point should advance";
			ASSERT_EQ(buffer.Remaining(), maxLength - readLength) << "Remaining should be decreased";
			ASSERT_THROW(buffer.ReadBits(-1), std::range_error*);
			ASSERT_THROW(buffer.ReadBits(65), std::range_error*);
			if ((value & 1ULL) == 0) {
				value |= 1ULL;
			} else {
				value = (value << 1) ^ value;
			}
		}
	}
	ASSERT_EQ(buffer.ReadByte(), static_cast<uint8_t>(0)) << "Reads beyond end should return 0";
	ASSERT_EQ(buffer.ReadShort(), static_cast<uint16_t>(0)) << "Reads beyond end should return 0";
	ASSERT_EQ(buffer.ReadInt(), static_cast<uint32_t>(0)) << "Reads beyond end should return 0";
	ASSERT_EQ(buffer.ReadLong(), static_cast<uint64_t>(0)) << "Reads beyond end should return 0";
	size_t readLength2 = 0ULL;
	for (int width = 1; width <= 64; ++width) {
		uint64_t value = 0ULL;
		while (value < (1ULL << width)) {
			uint64_t readValue = buffer2.ReadBits(width);
			ASSERT_EQ(readValue, value) << "Readback should be of expected values";
			ASSERT_EQ(buffer2.Size(), maxLength) << "Size should not change";
			readLength2 += width;
			ASSERT_EQ(buffer2.ReadPos(), readLength2) << "Read point should advance";
			ASSERT_EQ(buffer2.Remaining(), maxLength - readLength2) << "Remaining should be decreased";
			ASSERT_THROW(buffer2.ReadBits(-1), std::range_error*);
			ASSERT_THROW(buffer2.ReadBits(65), std::range_error*);
			if ((value & 1ULL) == 0) {
				value |= 1ULL;
			} else {
				value = (value << 1) ^ value;
			}
		}
	}
}

TEST(BitBufferTests, OtherTests) {
	bitbuffer::BitBuffer buffer;
	buffer.WriteInt(100);
	buffer.WriteInt(200);
	buffer.WriteInt(300);
	buffer.WriteBits(1,1);
	bitbuffer::BitBuffer buffer2;
	buffer2.WriteInt(400);
	buffer2.WriteInt(500);
	buffer2.WriteInt(600);
	buffer2.WriteBits(3, 2);
	size_t prevSize = buffer.Size();
	buffer.Append(buffer2);
	ASSERT_EQ(buffer.Size(), prevSize + buffer2.Size()) << "size is concatenated";
	ASSERT_EQ(buffer.ReadInt(), 100) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 200) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 300) << "values match";
	ASSERT_EQ(buffer.ReadBits(1), 1ULL) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 400) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 500) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 600) << "values match";
	ASSERT_EQ(buffer.ReadBits(2), 3ULL) << "values match";
	prevSize = buffer.Size();
	buffer.ResetReadPtr();
	ASSERT_EQ(buffer.ReadPos(), 0ULL) << "reset moves pointer back";
	ASSERT_EQ(buffer.Size(), prevSize) << "reset doesn't change size";
	ASSERT_EQ(buffer.ReadInt(), 100) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 200) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 300) << "values match";
	buffer.SkipBits(1); //exercise SkipBits
	ASSERT_EQ(buffer.ReadInt(), 400) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 500) << "values match";
	ASSERT_EQ(buffer.ReadInt(), 600) << "values match";
	ASSERT_EQ(buffer.ReadBits(2), 3ULL) << "values match";
	buffer.Clear();
	ASSERT_EQ(buffer.Size(), 0ULL);
	ASSERT_EQ(buffer.ReadPos(), 0ULL);
	ASSERT_EQ(buffer.Remaining(), 0ULL);
	ASSERT_EQ(buffer.ReadInt(), 0ULL);
}

TEST(BitBufferTests, ZigZagTests) {
	for (int32_t x = -80000; x < 80000 / 2; ++x) {
		uint32_t enc = bitbuffer::zigzagEncode(x);
		ASSERT_LE(std::bit_width(enc), 2 * std::bit_width(static_cast<uint32_t>(x))) << "bit size is bounded";
		int32_t dec = bitbuffer::zigzagDecode(enc);
		ASSERT_EQ(dec, x) << "encoding is reversible";
	}
	for (int32_t x = std::numeric_limits<int32_t>::min() / 2; x < std::numeric_limits<int32_t>::max() / 2; x += 0x10000) {
		uint32_t enc = bitbuffer::zigzagEncode(x);
		ASSERT_LE(std::bit_width(enc), 2 * std::bit_width(static_cast<uint32_t>(x))) << "bit size is bounded";
		int32_t dec = bitbuffer::zigzagDecode(enc);
		ASSERT_EQ(dec, x) << "encoding is reversible";
	}
}

TEST(BitBufferTests, GolombTests) {
	for (uint32_t m = 1; m < 256; ++m) {
		bitbuffer::BitBuffer buffer;
		for (uint32_t x = 0; x < std::numeric_limits<uint16_t>::max(); x += (x < 3*m) ? 1 : 997) {
			size_t prevSize = buffer.Size();
			bitbuffer::writeGolombCode(x, m, buffer);
			size_t newBits = buffer.Size() - prevSize;
			ASSERT_EQ(bitbuffer::golombCodeLength(x, m), newBits) << "length estimate is correct " << x << ", " << m;
			uint32_t dec = bitbuffer::readGolombCode(m, buffer);
			ASSERT_EQ(dec, x) << "Readback is OK";
			ASSERT_EQ(buffer.Remaining(), 0ULL) << "readback fully consumes written";
		}
	}
}

TEST(BitBufferTests, EliasTests) {
	bitbuffer::BitBuffer buffer;
	for (uint32_t x = 0; x < std::numeric_limits<uint16_t>::max();  x += 997) {
		size_t prevSize = buffer.Size();
		bitbuffer::writeEliasCode(x, buffer);
		size_t newBits = buffer.Size() - prevSize;
		ASSERT_EQ(bitbuffer::eliasCodeLength(x), newBits) << "length estimate is correct";
		uint32_t dec = bitbuffer::readEliasCode(buffer);
		ASSERT_EQ(dec, x) << "Readback is OK";
		ASSERT_EQ(buffer.Remaining(), 0ULL) << "readback fully consumes written";
	}
}

TEST(BitBufferTests, EliasFanoSequenceTest) {
	bitbuffer::BitBuffer buffer;
	std::vector<uint16_t> test = {
		1, 1, 2, 4, 6, 7, 8, 8, 9, 10, 13, 15, 15
	};
	bitbuffer::writeEliasFanoSequenceCode(test, 15, buffer);
	EXPECT_GE(bitbuffer::eliasFanoSequenceCodeLength(test.size(), 15), static_cast<uint32_t>(buffer.Size()));
	std::vector<uint16_t> readback = bitbuffer::readEliasFanoSequenceCode(test.size(), 15, buffer);
	ASSERT_EQ(readback, test);
	buffer.Clear();
	bitbuffer::writeEliasFanoSequenceCode(test, 255, buffer);
	EXPECT_GE(bitbuffer::eliasFanoSequenceCodeLength(test.size(), 255), static_cast<uint32_t>(buffer.Size()));
	std::vector<uint16_t> readback2 = bitbuffer::readEliasFanoSequenceCode(test.size(), 255, buffer);
	ASSERT_EQ(readback2, test);
}

TEST(BitBufferTests, EliasFanoSequenceTest2) {
	bitbuffer::BitBuffer buffer;
	std::vector<uint16_t> test;
	std::mt19937 rand;
	uint16_t i = 0;
	const uint16_t maxValToTest = 1000;
	while (i < maxValToTest) {
		i += rand() % 3;
		if (i < maxValToTest) {
			test.push_back(i);
		}
	}
	bitbuffer::writeEliasFanoSequenceCode(test, maxValToTest, buffer);
	EXPECT_GE(bitbuffer::eliasFanoSequenceCodeLength(test.size(), maxValToTest), static_cast<uint32_t>(buffer.Size()));
	std::vector<uint16_t> readback = bitbuffer::readEliasFanoSequenceCode(test.size(), maxValToTest, buffer);
	ASSERT_EQ(readback, test);
}

TEST(BitBufferTests, EliasFanoSequenceTest3) {
	bitbuffer::BitBuffer buffer;
	std::vector<uint16_t> test;
	std::mt19937 rand;
	uint16_t i = 0;
	const uint16_t maxValToTest = 1000;
	while (i < maxValToTest) {
		i += 5 + rand() % 5;
		if (i < maxValToTest) {
			test.push_back(i);
		}
	}
	bitbuffer::writeEliasFanoSequenceCode(test, maxValToTest, buffer);
	EXPECT_GE(bitbuffer::eliasFanoSequenceCodeLength(test.size(), maxValToTest), static_cast<uint32_t>(buffer.Size()));
	std::vector<uint16_t> readback = bitbuffer::readEliasFanoSequenceCode(test.size(), maxValToTest, buffer);
	ASSERT_EQ(readback, test);
}