#include "pch.h"

#include "../CompressionLib/inc/Huffman.h"
#include <random>

TEST(HuffmanTests, EmptyTest) {
	bitbuffer::BitBuffer buffer;
	std::vector<uint16_t> empty;
	huffman::huffmanEncode(empty, buffer);
	std::vector<uint16_t> readback;
	huffman::huffmanDecode(buffer, readback);
	ASSERT_EQ(readback, empty);
	ASSERT_EQ(buffer.Remaining(), 0ULL); //consumes all bits
}

TEST(HuffmanTests, BasicTest) {
	// use algorithm to generate all integer partitions to use as symbol frequencies from https://jeromekelleher.net/category/combinatorics.html
	// which should ensure lots of different huffman tree shapes
	const size_t n = 10;
	std::mt19937 rand;
	std::vector<uint16_t> counts(n + 1);
	size_t k = 1;
	counts[0] = 0;
	counts[1] = n;
	while (k != 0) {
		uint16_t x = counts[k - 1] + 1;
		uint16_t y = counts[k] - 1;
		--k;
		while (x <= y) {
			counts[k] = x;
			y -= x;
			++k;
		}
		counts[k] = x + y;
		std::vector<uint16_t> test;
		for (int i = 0; i <= k; ++i) {
			for (uint16_t count = 0; count < counts[i]; ++count) {
				test.push_back(i);
			}
		}
		std::shuffle(test.begin(), test.end(), rand);
		bitbuffer::BitBuffer buffer;
		huffman::huffmanEncode(test, buffer);
		std::vector<uint16_t> readback;
		huffman::huffmanDecode(buffer, readback);
		ASSERT_EQ(readback, test);
		ASSERT_EQ(buffer.Remaining(), 0ULL); //consumes all bits
	}
}

TEST(HuffmanTests, BasicTest2) {
	bitbuffer::BitBuffer buffer;
	std::vector<uint16_t> test = {
		0,
		0,
		0,
		1,
		1,
		2,
		2,
		3,
		3,
		4,
		4,
		15 // low probability all-bits value to test for clash with logic of pseudo-eof
	};
	huffman::huffmanEncode(test, buffer);
	std::string data = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Eleifend quam adipiscing vitae proin sagittis nisl rhoncus mattis rhoncus. Tellus orci ac auctor augue mauris augue neque gravida in. Elit scelerisque mauris pellentesque pulvinar pellentesque habitant morbi. Et malesuada fames ac turpis. Ullamcorper a lacus vestibulum sed arcu non. Tortor id aliquet lectus proin. Ut porttitor leo a diam sollicitudin tempor id eu.";
	std::vector<uint16_t> test2(data.cbegin(), data.cend());
	huffman::huffmanEncode(test2, buffer);
	std::vector<uint16_t> readback;
	huffman::huffmanDecode(buffer, readback);
	ASSERT_EQ(readback, test);
	std::vector<uint16_t> readback2;
	huffman::huffmanDecode(buffer, readback2);
	ASSERT_EQ(readback2, test2);
	ASSERT_EQ(buffer.Remaining(), 0ULL); //consumes all bits
}


TEST(HuffmanTests, BasicTest3) {
	bitbuffer::BitBuffer buffer;
	std::vector<uint16_t> test;
	buffer.WriteBits(0x1555, 13); //put some arbitrary data into buffer
	huffman::huffmanEncode(test, buffer);
	std::vector<uint16_t> readback;
	buffer.SkipBits(13);
	huffman::huffmanDecode(buffer, readback);
	ASSERT_EQ(readback, test);
	ASSERT_EQ(buffer.Remaining(), 0ULL); //consumes all bits
}

TEST(HuffmanTests, LargeTest) {
	std::mt19937 rand;
	bitbuffer::BitBuffer buffer;
	const size_t n = 100000ULL;
	std::vector<uint16_t> test;
	test.reserve(n);
	for (size_t i = 0ULL; i < n; ++i) {
		if ((rand() & 1) == 0) {
			test.push_back(rand() % 0xFFFF);
		} else {
			test.push_back(100 + (rand() % 5)); // make 100-105 more common
		}
	}
	huffman::huffmanEncode(test, buffer);
	std::vector<uint16_t> readback;
	huffman::huffmanDecode(buffer, readback);
	ASSERT_EQ(readback, test);
	ASSERT_EQ(buffer.Remaining(), 0ULL); //consumes all bits
}

TEST(HuffmanTests, CorruptStreamTest) {
	bitbuffer::BitBuffer buffer;
	std::vector<uint16_t> out;
	ASSERT_THROW(huffman::huffmanDecode(buffer, out), std::range_error*);
	buffer.WriteBits(1,1);
	ASSERT_THROW(huffman::huffmanDecode(buffer, out), std::range_error*);
	buffer.Clear();
	std::string data = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Eleifend quam adipiscing vitae proin sagittis nisl rhoncus mattis rhoncus. Tellus orci ac auctor augue mauris augue neque gravida in. Elit scelerisque mauris pellentesque pulvinar pellentesque habitant morbi. Et malesuada fames ac turpis. Ullamcorper a lacus vestibulum sed arcu non. Tortor id aliquet lectus proin. Ut porttitor leo a diam sollicitudin tempor id eu.";
	std::vector<uint16_t> test(data.cbegin(), data.cend());
	huffman::huffmanEncode(test, buffer);
	buffer.SkipBits(1);
	ASSERT_THROW(huffman::huffmanDecode(buffer, out), std::range_error*);
}