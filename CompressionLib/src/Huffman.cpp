#include "../inc/Huffman.h"
#include <iostream>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <format>

namespace huffman {
	struct HuffmanEntry {
		uint32_t symbol;
		size_t frequency;
		uint8_t codeBitLength;
		uint32_t code;
	};

	struct HuffmanNode {
		uint8_t depth;
		size_t frequency;
		std::vector<uint32_t> children;
	};

	struct NodeCompare {
		bool operator()(const HuffmanNode& a, const HuffmanNode& b) {
			return a.frequency > b.frequency;
		}
	};

	struct CanonicalSorter {
		bool operator()(const HuffmanEntry& a, const HuffmanEntry& b) {
			if (a.codeBitLength == b.codeBitLength) {
				return a.symbol < b.symbol;
			} else {
				return a.codeBitLength < b.codeBitLength;
			}
		}
	};

	const uint32_t PSEUDO_EOF = 0xFFFFFFFFU;

	std::unique_ptr<uint8_t[]> huffmanEncode(const std::vector<uint16_t>& data, size_t& outputLength) {
		bitbuffer::BitBuffer buffer;
		huffmanEncode(data, buffer);
		return buffer.Save(outputLength);
	}

	void huffmanEncode(const std::vector<uint16_t>& data, bitbuffer::BitBuffer& buffer) {
		std::unordered_map<uint32_t, HuffmanEntry> freqs;
		uint8_t symbolBitWidth = 1;
		for (uint16_t symbol : data) {
			uint32_t intSymbol = static_cast<uint32_t>(symbol);
			uint8_t bitWidth = std::bit_width(intSymbol);
			if (bitWidth > symbolBitWidth) {
				symbolBitWidth = bitWidth;
			}
			if (!freqs.contains(intSymbol)) {
				freqs[symbol] = HuffmanEntry{
					.symbol = intSymbol,
					.frequency = 0ULL,
					.codeBitLength = 0,
					.code = 0
				};
			}
			++freqs[intSymbol].frequency;
		}
		//add pseudo-EOF
		freqs[PSEUDO_EOF] = HuffmanEntry{
					.symbol = PSEUDO_EOF,
					.frequency = 0ULL,
					.codeBitLength = static_cast<uint8_t>(data.empty() ? 1 : 0),
					.code = 0
		};
		std::priority_queue<HuffmanNode, std::vector<HuffmanNode>, NodeCompare> depthBuilder;
		for (const auto& leaf : freqs) {
			std::vector<uint32_t> child;
			child.push_back(leaf.first);
			depthBuilder.push(HuffmanNode{
				.depth = 0,
				.frequency = leaf.second.frequency,
				.children = child
				});
		}
		while (depthBuilder.size() > 1) {
			const HuffmanNode& a = depthBuilder.top();
			size_t totalFreq = a.frequency;
			uint8_t maxDepth = a.depth;
			std::vector<uint32_t> combinedChildren(a.children.cbegin(), a.children.cend());
			depthBuilder.pop();
			const HuffmanNode& b = depthBuilder.top();
			totalFreq += b.frequency;
			maxDepth = std::max(maxDepth, b.depth);
			combinedChildren.insert(combinedChildren.end(), b.children.cbegin(), b.children.cend());
			depthBuilder.pop();
			for (uint32_t child : combinedChildren) {
				++freqs[child].codeBitLength;
			}
			++maxDepth;
			HuffmanNode newNode{
				.depth = maxDepth,
				.frequency = totalFreq,
				.children = std::move(combinedChildren)
			};
			depthBuilder.push(newNode);
		}
		std::vector<HuffmanEntry> sortedEntries;
		sortedEntries.reserve(freqs.size());
		std::transform(freqs.begin(), freqs.end(), std::back_inserter(sortedEntries), [](auto& kv) { return kv.second; });
		std::sort(sortedEntries.begin(), sortedEntries.end(), CanonicalSorter());
		std::unordered_map<uint32_t, HuffmanEntry> symbolMap;
		uint8_t maxLength = std::max(depthBuilder.top().depth, static_cast<uint8_t>(1));
		buffer.WriteByte(maxLength);
		std::vector<uint16_t> lengths;
		uint8_t prevLength = 0;
		uint32_t code = 0;
		uint32_t count = 0;
		for (HuffmanEntry& entry : sortedEntries) {
			if (entry.codeBitLength != prevLength) {
				if (prevLength != 0) {
					buffer.WriteShort(count);
					lengths.push_back(count);
				}
				for (uint8_t length = prevLength + 1; length < entry.codeBitLength; ++length) {
					buffer.WriteShort(0);
					lengths.push_back(0);
				}
				count = 0;
				code = code << (entry.codeBitLength - prevLength);
			}
			entry.code = code;
			++code;
			++count;
			prevLength = entry.codeBitLength;
			symbolMap[entry.symbol] = entry;
		}
		buffer.WriteShort(count);
		lengths.push_back(count);
		buffer.WriteByte(symbolBitWidth);
		size_t lengthIndex = 0;
		uint16_t mask = ((1U << symbolBitWidth) - 1U);
		std::vector<HuffmanEntry>::const_iterator symbolIter = sortedEntries.cbegin();
		while (symbolIter != sortedEntries.cend()) {
			uint16_t rangeLength = lengths[lengthIndex++];
			std::vector<uint16_t> symbolRange;
			symbolRange.reserve(rangeLength);
			std::transform(symbolIter,
				symbolIter + rangeLength,
				std::back_inserter(symbolRange),
				[&](auto& entry) { return static_cast<uint16_t>(entry.symbol) & mask; });
			if (bitbuffer::eliasFanoSequenceCodeLength(rangeLength, mask) < static_cast<uint32_t>(rangeLength) * static_cast<uint32_t>(symbolBitWidth)) {
				bitbuffer::writeEliasFanoSequenceCode(symbolRange, mask, buffer);
			} else {
				for (const uint16_t& sym : symbolRange) {
					buffer.WriteBits(sym, symbolBitWidth);
				}
			}
			symbolIter += rangeLength;
		}
		for (uint16_t symbol : data) {
			uint32_t intSymbol = static_cast<uint32_t>(symbol);
			const HuffmanEntry& entry = symbolMap[intSymbol];
			buffer.WriteBits(entry.code, entry.codeBitLength);
		}
		buffer.WriteBits(symbolMap[PSEUDO_EOF].code, symbolMap[PSEUDO_EOF].codeBitLength);
	}

	std::vector<uint16_t> huffmanDecode(const uint8_t encoded[], size_t compressedLength) {
		bitbuffer::BitBuffer buffer;
		buffer.Load(encoded, 0, 8ULL * compressedLength);
		std::vector<uint16_t> results;
		huffmanDecode(buffer, results);
		return results;
	}

	void huffmanDecode(bitbuffer::BitBuffer& buffer, std::vector<uint16_t>& decoded) {
		uint8_t maxLength = buffer.ReadByte();
		std::vector<uint16_t> lengthCounts;
		uint16_t totalSymbols = 0;
		for (uint8_t length = 1; length <= maxLength; ++length) {
			uint16_t count = buffer.ReadShort();
			totalSymbols += count;
			lengthCounts.push_back(count);
		}
		uint8_t symbolBitWidth = static_cast<uint8_t>(buffer.ReadByte());
		uint16_t mask = ((1U << symbolBitWidth) - 1U);
		std::vector<uint16_t> symbolList;
		symbolList.reserve(totalSymbols);
		for (uint16_t rangeLength : lengthCounts) {
			if (bitbuffer::eliasFanoSequenceCodeLength(rangeLength, mask) < static_cast<uint32_t>(rangeLength) * static_cast<uint32_t>(symbolBitWidth)) {
				std::vector<uint16_t> symbolRange = bitbuffer::readEliasFanoSequenceCode(rangeLength, mask, buffer);
				symbolList.insert(symbolList.end(), symbolRange.cbegin(), symbolRange.cend());
			} else {
				for (int i = 0; i < rangeLength; ++i) {
					symbolList.push_back(static_cast<uint16_t>(buffer.ReadBits(symbolBitWidth)));
				}
			}
		}
		std::unordered_map<uint32_t, HuffmanEntry> codeMap;
		uint8_t length = 1;
		uint32_t code = 0;
		uint8_t prevLength = 0;
		uint16_t symbolOffset = 0;
		for (uint16_t symbolCount = 0; symbolCount < totalSymbols; ++symbolCount) {
			while (lengthCounts[static_cast<size_t>(length) - 1] == 0) {
				++length;
			}
			--lengthCounts[static_cast<size_t>(length) - 1];
			uint32_t symbol = static_cast<uint32_t>(symbolList[symbolOffset++]);
			if (length != prevLength) {
				code = code << (length - prevLength);
			}
			if (symbolCount == totalSymbols - 1) {
				if (symbol != ((1U << symbolBitWidth) - 1U)) {
					throw new std::range_error("Invalid bitstream");
				}
				symbol = PSEUDO_EOF;
			}
			codeMap[code] = HuffmanEntry {
				.symbol = symbol,
				.frequency = 0,
				.codeBitLength = length,
				.code = code
			};
			++code;
			prevLength = length;
		}
		uint32_t decode = 0;
		uint8_t bitcount = 0;
		while(buffer.Remaining() > 0ULL) {
			decode = (decode << 1) | static_cast<uint32_t>(buffer.ReadBits(1));
			++bitcount;
			if (codeMap.contains(decode)) {
				const HuffmanEntry& entry = codeMap[decode];
				if (entry.codeBitLength == bitcount) {
					if (entry.symbol == PSEUDO_EOF) {
						return;
					}
					decoded.push_back(static_cast<uint16_t>(entry.symbol));
					decode = 0;
					bitcount = 0;
				}
			}
		}
		throw new std::range_error("Invalid bitstream");
	}


	std::vector<uint16_t> runLengthEncode(const std::vector<uint16_t>& data) {
		std::vector<uint16_t> retval;
		if (data.empty()) {
			return retval;
		}
		retval.reserve(data.size());
		uint16_t prevVal = 0;
		bool newRun = true;
		uint16_t count = 0;
		for (const uint16_t& val : data) {
			if (val == prevVal && !newRun) {
				++count;
				if (count == 1) {
					retval.push_back(val);
				} else if (count >= 0x8000) {
					retval.push_back(count - 1);
					count = 0;
					newRun = true;
				}
			} else {
				newRun = false;
				if (count > 0) {
					retval.push_back(count - 1);
					count = 0;
				}
				prevVal = val;
				retval.push_back(val);
			}
		}
		if (count > 0) {
			retval.push_back(count - 1);
		}
		return retval;
	}

	std::vector<uint16_t> runLengthDecode(const std::vector<uint16_t>& data) {
		std::vector<uint16_t> retval;
		if (data.empty()) {
			return retval;
		}
		retval.reserve(data.size());
		uint16_t prevVal = 0;
		bool pendingCount = false;
		bool newRun = true;
		for (const uint16_t& val : data) {
			if (pendingCount) {
				for (uint16_t i = 0; i < val; ++i) {
					retval.push_back(prevVal);
				}
				pendingCount = false;
				newRun = true;
			} else {
				retval.push_back(val);
				if (val == prevVal && !newRun) {
					pendingCount = true;
				}
				newRun = false;
				prevVal = val;
			}
		}
		return retval;
	}
}