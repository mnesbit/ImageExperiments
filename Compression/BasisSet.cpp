#include "BasisSet.h"
#include "../SimpleMatrix/inc/symmeigen.h"
#include <set>
#include <algorithm>

namespace basis {

	// fully substituted version of optimised summedLaplaceModel from FitCovarianceModel experiments
	double covariancePredictionModelY(double dx, double dy) {
		return 3817.7299999999996 * exp(-1.48854e-05 * dx * dx + -1.7273e-05 * dy * dy)
			+ 657.8100000000001 * exp(-0.0436057 * abs(dx) + -0.050844400000000005 * abs(dy));
	}
	// fully substituted version of optimised ellipticLaplaceModel for U channel
	double covariancePredictionModelU(double dx, double dy) {
		return 241.49 * exp(-0.00134755 * abs(dx) + -0.00147572 * abs(dy));
	}

	// fully substituted version of optimised ellipticLaplaceModel for V channel
	double covariancePredictionModelV(double dx, double dy) {
		return 371.87199999999996 * exp(-0.00147084 * abs(dx) + -0.0015265799999999998 * abs(dy));
	}

	std::vector<math::Vector> createBasis(math::Matrix& covariance)
	{
		math::SymmetricEigenDecomposition eigenDeco(covariance);
		math::Matrix eigenVectors = eigenDeco.EigenVectors();
		math::Vector eigenValues = eigenDeco.EigenValues();
		std::vector<std::pair<double, int> > sortedIndex;
		sortedIndex.reserve(eigenValues.Length());
		for (int i = 0; i < eigenValues.Length(); ++i) {
			sortedIndex.push_back(std::make_pair(abs(eigenValues[i]), i));
		}
		std::sort(sortedIndex.begin(), sortedIndex.end(), std::greater<std::pair<double, int> >());
		std::vector<math::Vector> basisSet;
		for (int i = 0; i < eigenValues.Length(); ++i) {
			math::Vector col = eigenVectors.GetColumn(sortedIndex[i].second);
			basisSet.push_back(col);
		}
		return basisSet;
	}

	typedef struct Pt_Tag
	{
		int x, y;
	} Pt;

	static double Dist(Pt& line1, Pt& line2, Pt& point) {
		return static_cast<double>((line2.x - line1.x) * (line1.y - point.y) - (line1.x - point.x) * (line2.y - line1.y)) / sqrt((line2.x - line1.x) * (line2.x - line1.x) + (line2.y - line1.y) * (line2.y - line1.y));
	}

	struct ShapeComparator {
		bool operator()(std::vector<double> a, std::vector<double> b) const
		{
			for (int i = 0; i < a.size(); ++i) {
				if (signbit(a[i]) != signbit(b[i])) {
					return signbit(b[i]);
				}
			}
			return false;
		}
	};

	std::vector<std::vector<double> > distinctLineShapes(size_t blockSize) {
		std::set<std::vector<double>, ShapeComparator> basisSet;
		basisSet.insert(std::vector<double>(blockSize * blockSize, true)); //dc basis
		for (int side1 = -static_cast<int>(blockSize); side1 < 2 * static_cast<int>(blockSize); ++side1)
		{
			Pt s1 = { .x = side1, .y = -static_cast<int>(blockSize) };
			for (int side2 = -static_cast<int>(blockSize); side2 < 2 * static_cast<int>(blockSize); ++side2)
			{
				Pt s2a = { .x = -static_cast<int>(blockSize), .y = side2 };
				Pt s2b = { .x = side2, .y = 2 * static_cast<int>(blockSize) };
				Pt s2c = { .x = 2 * static_cast<int>(blockSize), .y = side2 };
				std::vector<double> bits1(blockSize * blockSize);
				std::vector<double> invbits1(blockSize * blockSize);
				std::vector<double> bits2(blockSize * blockSize);
				std::vector<double> invbits2(blockSize * blockSize);
				std::vector<double> bits3(blockSize * blockSize);
				std::vector<double> invbits3(blockSize * blockSize);
				for (size_t x = 0; x < blockSize; ++x)
				{
					for (size_t y = 0; y < blockSize; ++y)
					{
						Pt s3 = { .x = static_cast<int>(x), .y = static_cast<int>(y) };
						bits1[x + y * blockSize] = tanh(2.0 * Dist(s1, s2a, s3));
						invbits1[x + y * blockSize] = -bits1[x + y * blockSize];
						bits2[x + y * blockSize] = tanh(2.0 * Dist(s1, s2b, s3));
						invbits2[x + y * blockSize] = -bits2[x + y * blockSize];
						bits3[x + y * blockSize] = tanh(2.0 * Dist(s1, s2c, s3));
						invbits3[x + y * blockSize] = -bits3[x + y * blockSize];
					}
				}
				if (!(basisSet.contains(bits1) || basisSet.contains(invbits1)))
				{
					basisSet.insert(bits1);
				}
				if (!(basisSet.contains(bits2) || basisSet.contains(invbits2)))
				{
					basisSet.insert(bits2);
				}
				if (!(basisSet.contains(bits3) || basisSet.contains(invbits3)))
				{
					basisSet.insert(bits3);
				}
			}
		}
		for (int side1 = -static_cast<int>(blockSize); side1 < 2 * static_cast<int>(blockSize); ++side1)
		{
			Pt s1 = { .x = 2 * static_cast<int>(blockSize) , .y = side1 };
			for (int side2 = -static_cast<int>(blockSize); side2 < 2 * static_cast<int>(blockSize); ++side2)
			{
				Pt s2a = { .x = -static_cast<int>(blockSize), .y = side2 };
				Pt s2b = { .x = side2, .y = 2 * static_cast<int>(blockSize) };
				Pt s2c = { .x = -static_cast<int>(blockSize), .y = side1 };
				std::vector<double> bits1(blockSize * blockSize);
				std::vector<double> invbits1(blockSize * blockSize);
				std::vector<double> bits2(blockSize * blockSize);
				std::vector<double> invbits2(blockSize * blockSize);
				std::vector<double> bits3(blockSize * blockSize);
				std::vector<double> invbits3(blockSize * blockSize);
				for (size_t x = 0; x < blockSize; x++)
				{
					for (size_t y = 0; y < blockSize; y++)
					{
						Pt s3{ .x = static_cast<int>(x), .y = static_cast<int>(y) };
						bits1[x + y * blockSize] = tanh(2.0 * Dist(s1, s2a, s3));
						invbits1[x + y * blockSize] = -bits1[x + y * blockSize];
						bits2[x + y * blockSize] = tanh(2.0 * Dist(s1, s2b, s3));
						invbits2[x + y * blockSize] = -bits2[x + y * blockSize];
						bits3[x + y * blockSize] = tanh(2.0 * Dist(s2b, s2c, s3));
						invbits3[x + y * blockSize] = -bits3[x + y * blockSize];
					}
				}
				if (!(basisSet.contains(bits1) || basisSet.contains(invbits1)))
				{
					basisSet.insert(bits1);
				}
				if (!(basisSet.contains(bits2) || basisSet.contains(invbits2)))
				{
					basisSet.insert(bits2);
				}
				if (!(basisSet.contains(bits3) || basisSet.contains(invbits3)))
				{
					basisSet.insert(bits3);
				}
			}
		}
		return std::vector(basisSet.cbegin(), basisSet.cend());
	}

	math::Matrix createSegmentDictionary(size_t blockSize, const std::vector<std::vector<double> >& basisSet) {
		math::Matrix dictionary(basisSet.size(), blockSize * blockSize);
		size_t i = 0;
		for (const std::vector<double>& bits : basisSet) {
			double total = 0.0;
			bool allSet = true;
			bool allClear = true;
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				if (bits[j] >= 0.0) {
					allClear = false;
				} else {
					allSet = false;
				}
				dictionary[i][j] = bits[j];
				total += bits[j];
			}
			double mean = total / static_cast<double>(blockSize * blockSize);
			double sumsq = 0.0;
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				double value = dictionary[i][j];
				if (!(allSet || allClear)) {
					value -= mean;
				}
				dictionary[i][j] = value;
				sumsq += value * value;
			}
			double norm = sqrt(sumsq);
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				double value = dictionary[i][j];
				if (sumsq != 0.0) {
					dictionary[i][j] = value / norm;
				}
			}
			++i;
		}
		return dictionary;
	}

	math::Matrix createIntraSegmentDictionary(size_t blockSize, const std::vector<double>& segmentMask, CovarianceModel model) {
		int count1 = 0;
		int count2 = 0;
		math::Matrix part1(blockSize * blockSize, blockSize * blockSize);
		math::Matrix part2(blockSize * blockSize, blockSize * blockSize);
		for (size_t x = 0; x < blockSize * blockSize; ++x) {
			int x1 = static_cast<int>(x % blockSize);
			int y1 = static_cast<int>(x / blockSize);
			if (segmentMask[x] >= 0.0) ++count1;
			if (segmentMask[x] <= 0.0) ++count2;
			for (size_t y = 0; y < blockSize * blockSize; ++y) {
				int x2 = static_cast<int>(y % blockSize);
				int y2 = static_cast<int>(y / blockSize);
				double distx = static_cast<double>(x2 - x1);
				double disty = static_cast<double>(y2 - y1);
				part1[x][y] = std::max(0.0, 0.5 + segmentMask[x]) * std::max(0.0, 0.5 + segmentMask[y]) * model(distx, disty);
				part2[x][y] = std::max(0.0, 0.5 - segmentMask[x]) * std::max(0.0, 0.5 - segmentMask[y]) * model(distx, disty);
			}
		}

		std::vector<math::Vector> basis1 = createBasis(part1);
		std::vector<math::Vector> basis2 = createBasis(part2);
		size_t part1Count = std::max(0, count1 - 2);
		size_t part2Count = std::max(0, count2 - 2);

		math::Matrix dictionary(part1Count + part2Count, blockSize * blockSize);
		for (size_t i = 0; i < part1Count; ++i) {
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				dictionary[i][j] = basis1[i + 1][j];
			}
		}
		for (size_t i = 0; i < part2Count; ++i) {
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				dictionary[part1Count + i][j] = basis2[i + 1][j];
			}
		}
		return dictionary;
	}
}