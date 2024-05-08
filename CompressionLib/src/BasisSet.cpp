#include "../inc/BasisSet.h"
#include "../../SimpleMatrix/inc/symmeigen.h"
#include <map>
#include <numbers>
#include <algorithm>
#include <numbers>
#include "../../eigen/Eigen/Eigenvalues"

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

	math::Matrix createDCTDictionary(int blockSize) {
		math::Matrix dctDict(blockSize * blockSize, blockSize * blockSize);
		for (int u0 = 0; u0 < blockSize; ++u0) {
			for (int v0 = 0; v0 < blockSize; ++v0) {
				int index = (u0 * blockSize) + v0;
				for (int x = 0; x < blockSize; ++x) {
					for (int y = 0; y < blockSize; ++y) {
						double sum = cos(((2 * x + 1) / (2.0 * blockSize)) * u0 * std::numbers::pi) * cos(((2 * y + 1) / (2.0 * blockSize)) * v0 * std::numbers::pi);
						if (u0 && v0) {
							sum = sum * sqrt(4.0 / (blockSize * blockSize));
						} else if (u0 || v0) {
							sum = sum * sqrt(2.0 / (blockSize * blockSize));
						} else {
							sum = 0.5 * sum * sqrt(4.0 / (blockSize * blockSize));
						}
						dctDict[index][x + (blockSize * y)] = sum;
					}
				}
			}
		}
		return dctDict;
	}

	Eigen::MatrixXf createDCTDictionaryFast(int blockSize) {
		Eigen::MatrixXf dctDict(blockSize * blockSize, blockSize * blockSize);
		for (int u0 = 0; u0 < blockSize; ++u0) {
			for (int v0 = 0; v0 < blockSize; ++v0) {
				int index = (u0 * blockSize) + v0;
				for (int x = 0; x < blockSize; ++x) {
					for (int y = 0; y < blockSize; ++y) {
						double sum = cos(((2 * x + 1) / (2.0 * blockSize)) * u0 * std::numbers::pi) * cos(((2 * y + 1) / (2.0 * blockSize)) * v0 * std::numbers::pi);
						if (u0 && v0) {
							sum = sum * sqrt(4.0 / (blockSize * blockSize));
						} else if (u0 || v0) {
							sum = sum * sqrt(2.0 / (blockSize * blockSize));
						} else {
							sum = 0.5 * sum * sqrt(4.0 / (blockSize * blockSize));
						}
						dctDict(index, x + (blockSize * y)) = static_cast<float>(sum);
					}
				}
			}
		}
		return dctDict;
	}

	math::Matrix createKLTDictionary(int blockSize, CovarianceModel model) {
		math::Matrix covarianceMatrix(blockSize * blockSize, blockSize * blockSize);
		for (int i = 0; i < blockSize * blockSize; ++i) {
			int xpos1 = i / blockSize;
			int ypos1 = i % blockSize;
			for (int j = 0; j < blockSize * blockSize; ++j) {
				int xpos2 = j / blockSize;
				int ypos2 = j % blockSize;
				covarianceMatrix[i][j] = model(xpos1 - xpos2, ypos1 - ypos2);
			}
		}
		return createDictionary(covarianceMatrix);
	}

	Eigen::MatrixXf createKLTDictionaryFast(int blockSize, CovarianceModel model) {
		Eigen::MatrixXf covarianceMatrix(blockSize * blockSize, blockSize * blockSize);
		for (int i = 0; i < blockSize * blockSize; ++i) {
			int xpos1 = i / blockSize;
			int ypos1 = i % blockSize;
			for (int j = 0; j < blockSize * blockSize; ++j) {
				int xpos2 = j / blockSize;
				int ypos2 = j % blockSize;
				covarianceMatrix(i,j) = static_cast<float>(model(xpos1 - xpos2, ypos1 - ypos2));
			}
		}
		return createDictionaryFast(covarianceMatrix);
	}

	math::Matrix createDictionary(const math::Matrix& covariance) {
		std::vector<math::Vector> basisVec = createBasis(covariance);
		math::Matrix dictionary(basisVec.size(), covariance.Columns());
		for (size_t i = 0; i < covariance.Columns(); ++i) {
			memcpy(dictionary.Data() + (i * covariance.Columns()), basisVec[i].Data(), sizeof(double) * basisVec[i].Length());
		}
		return dictionary;
	}

	Eigen::MatrixXf createDictionaryFast(const Eigen::MatrixXf& covariance) {
		std::vector<Eigen::VectorXf> basisVec = createBasisFast(covariance);
		Eigen::MatrixXf dictionary(basisVec.size(), covariance.cols());
		for (int i = 0; i < covariance.rows(); ++i) {
			dictionary.row(i) = basisVec[i];
		}
		return dictionary;
	}

	std::vector<math::Vector> createBasis(const math::Matrix& covariance)
	{
		if (covariance.Rows() != covariance.Columns()) {
			throw new std::range_error("Requires a square symmetric covariance matrix");
		}
		if (covariance.Rows() == 0ULL) {
			return std::vector<math::Vector>();
		}
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
			// scan for first significant entry and sign flip eigenvector if less than zero to ensure consistency of representation
			for (int j = 0; j < col.Length(); ++j) {
				if (abs(col[j]) > 1E-10) {
					if (col[j] < 0.0) {
						for (int k = 0; k < col.Length(); ++k) {
							col[k] = -col[k];
						}
					}
					break;
				}
			}
			basisSet.push_back(col);
		}
		return basisSet;
	}

	std::vector<Eigen::VectorXf> createBasisFast(const Eigen::MatrixXf& covariance)
	{
		if (covariance.rows() != covariance.cols()) {
			throw new std::range_error("Requires a square symmetric covariance matrix");
		}
		if (covariance.rows() == 0ULL) {
			return std::vector<Eigen::VectorXf>();
		}
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigenDeco(covariance.rows());
		eigenDeco.compute(covariance);
		std::vector<std::pair<double, int> > sortedIndex;
		sortedIndex.reserve(covariance.rows());
		for (int i = 0; i < covariance.rows(); ++i) {
			sortedIndex.push_back(std::make_pair(abs(eigenDeco.eigenvalues()(i)), i));
		}
		std::sort(sortedIndex.begin(), sortedIndex.end(), std::greater<std::pair<double, int> >());
		std::vector<Eigen::VectorXf> basisSet;
		for (int i = 0; i < sortedIndex.size(); ++i) {
			const Eigen::VectorXf& col = eigenDeco.eigenvectors().col(sortedIndex[i].second);
			// scan for first significant entry and sign flip eigenvector if less than zero to ensure consistency of representation
			for (int j = 0; j < col.size(); ++j) {
				if (abs(col[j]) > 1E-10) {
					if (col[j] < 0.0) {
						basisSet.push_back(-1.0 * col);
					} else {
						basisSet.push_back(col);
					}
					break;
				}
			};
		}
		return basisSet;
	}

	static double Dist(const Pt& line1, const Pt& line2, const Pt& point) {
		return static_cast<double>((line2.x - line1.x) * (line1.y - point.y) - (line1.x - point.x) * (line2.y - line1.y)) / sqrt((line2.x - line1.x) * (line2.x - line1.x) + (line2.y - line1.y) * (line2.y - line1.y));
	}

	struct ShapeComparator {
		bool operator()(std::vector<bool> a, std::vector<bool> b) const
		{
			for (int i = 0; i < a.size(); ++i) {
				if (a[i] != b[i]) {
					return a[i];
				}
			}
			return false;
		}
	};

	std::vector<Line> distinctLineShapes(size_t blockSize) {
		std::map<std::vector<bool>, Line, ShapeComparator> basisSet;
		// Add strict horizontal and vertical (plus implicit DC) as otherwise faintly sloped versions enter in general rotuine below
		for (int side = -1; side < static_cast<int>(blockSize)+1; ++side) {
			Line horiz = Line{ .a = Pt {.x = 0, .y = side }, .b = Pt {.x = static_cast<int>(blockSize), .y = side } };
			Line vert = Line{ .a = Pt {.x = side, .y = 0 }, .b = Pt {.x = side, .y =static_cast<int>(blockSize) } };
			std::vector<bool> horizBits(blockSize * blockSize);
			std::vector<bool> vertBits(blockSize * blockSize);
			for (size_t x = 0; x < blockSize; ++x) {
				for (size_t y = 0; y < blockSize; ++y) {
					Pt s = { .x = static_cast<int>(x), .y = static_cast<int>(y) };
					horizBits[x + y * blockSize] = (Dist(horiz.a, horiz.b, s) >= 0.0);
					vertBits[x + y * blockSize] = (Dist(vert.a, vert.b, s) >= 0.0);
				}
			}
			if (!(basisSet.contains(horizBits))) {
				basisSet[horizBits] = horiz;
			}
			if (!(basisSet.contains(vertBits))) {
				basisSet[vertBits] = vert;
			}
		}
		for (int side1 = -static_cast<int>(blockSize);  side1 < 2 * static_cast<int>(blockSize); ++side1) {
			Pt s1 = { .x = side1, .y = -static_cast<int>(blockSize) };
			for (int side2 = -static_cast<int>(blockSize); side2 < 2 * static_cast<int>(blockSize); ++side2) {
				Pt s2a = { .x = -static_cast<int>(blockSize), .y = side2 };
				Pt s2b = { .x = side2, .y = static_cast<int>(blockSize) };
				Pt s2c = { .x = static_cast<int>(blockSize), .y = side2 };
				std::vector<bool> bits1(blockSize * blockSize);
				std::vector<bool> invbits1(blockSize * blockSize);
				std::vector<bool> bits2(blockSize * blockSize);
				std::vector<bool> invbits2(blockSize * blockSize);
				std::vector<bool> bits3(blockSize * blockSize);
				std::vector<bool> invbits3(blockSize * blockSize);
				for (size_t x = 0; x < blockSize; ++x) {
					for (size_t y = 0; y < blockSize; ++y) {
						Pt s3 = { .x = static_cast<int>(x), .y = static_cast<int>(y) };
						bits1[x + y * blockSize] = (Dist(s1, s2a, s3) >= 0.0);
						invbits1[x + y * blockSize] = !bits1[x + y * blockSize];
						bits2[x + y * blockSize] = (Dist(s1, s2b, s3) >= 0.0);
						invbits2[x + y * blockSize] = !bits2[x + y * blockSize];
						bits3[x + y * blockSize] = (Dist(s1, s2c, s3) >= 0.0);
						invbits3[x + y * blockSize] = !bits3[x + y * blockSize];
					}
				}
				if (!(basisSet.contains(bits1) || basisSet.contains(invbits1))) {
					basisSet[bits1] = Line{ .a = s1, .b = s2a };
				}
				if (!(basisSet.contains(bits2) || basisSet.contains(invbits2))) {
					basisSet[bits2] = Line{ .a = s1, .b = s2b };
				}
				if (!(basisSet.contains(bits3) || basisSet.contains(invbits3))) {
					basisSet[bits3] = Line{ .a = s1, .b = s2c };
				}
			}
		}
		for (int side1 = -static_cast<int>(blockSize); side1 < 2 * static_cast<int>(blockSize); ++side1) {
			Pt s1a = { .x = static_cast<int>(blockSize) , .y = side1 };
			Pt s1b = { .x = -static_cast<int>(blockSize), .y = side1 };
			for (int side2 = -static_cast<int>(blockSize); side2 < 2 * static_cast<int>(blockSize); ++side2) {
				Pt s2a = { .x = -static_cast<int>(blockSize), .y = side2 };
				Pt s2b = { .x = side2, .y = static_cast<int>(blockSize) };
				std::vector<bool> bits1(blockSize * blockSize);
				std::vector<bool> invbits1(blockSize * blockSize);
				std::vector<bool> bits2(blockSize * blockSize);
				std::vector<bool> invbits2(blockSize * blockSize);
				std::vector<bool> bits3(blockSize * blockSize);
				std::vector<bool> invbits3(blockSize * blockSize);
				for (size_t x = 0; x < blockSize; x++) {
					for (size_t y = 0; y < blockSize; y++) {
						Pt s3{ .x = static_cast<int>(x), .y = static_cast<int>(y) };
						bits1[x + y * blockSize] = (Dist(s1a, s2a, s3) >= 0.0);
						invbits1[x + y * blockSize] = !bits1[x + y * blockSize];
						bits2[x + y * blockSize] = (Dist(s1a, s2b, s3) >= 0.0);
						invbits2[x + y * blockSize] = !bits2[x + y * blockSize];
						bits3[x + y * blockSize] = (Dist(s1a, s2b, s3) >= 0.0);
						invbits3[x + y * blockSize] = !bits3[x + y * blockSize];
					}
				}
				if (!(basisSet.contains(bits1) || basisSet.contains(invbits1))) {
					basisSet[bits1] = Line{ .a = s1a, .b = s2a };
				}
				if (!(basisSet.contains(bits2) || basisSet.contains(invbits2))) {
					basisSet[bits2] = Line{ .a = s1a, .b = s2b };
				}
				if (!(basisSet.contains(bits3) || basisSet.contains(invbits3))) {
					basisSet[bits3] = Line{ .a = s1a, .b = s2b };
				}
			}
		}
		std::vector<Line> retval;
		std::transform(basisSet.cbegin(), basisSet.cend(), std::back_inserter(retval), [](const auto& x) { return x.second; });
		return retval;
	}

	math::Matrix createSegmentDictionary(size_t blockSize, const std::vector<Line>& basisSet) {
		math::Matrix dictionary(basisSet.size(), blockSize * blockSize);
		const double sigma = 1.0;
		const int kernelWidth = static_cast<int>(1 + sigma * 6);
		const int kernelHalfWidth = kernelWidth / 2;
		math::Vector gaussian(kernelWidth * kernelWidth);
		for (int dx = -kernelHalfWidth; dx <= +kernelHalfWidth; ++dx) {
			for (int dy = -kernelHalfWidth; dy <= +kernelHalfWidth; ++dy) {
				gaussian[(dx + kernelHalfWidth) + kernelWidth * (dy + kernelHalfWidth)] = exp(-static_cast<double>(dx * dx + dy * dy) / (2.0 * sigma * sigma));
			}
		}
		size_t i = 0;
		for (const Line separator : basisSet) {
			bool allSet = true;
			bool allClear = true;
			math::Vector zoomIn(4 * blockSize * blockSize);
			Line doubledLine = Line{
				.a = Pt {
					.x = separator.a.x * 2,
					.y = separator.a.y * 2
				},
				.b = Pt {
					.x = separator.b.x * 2,
					.y = separator.b.y * 2
				}
			};
			for (int x = 0; x < 2 * blockSize; ++x) {
				for (int y = 0; y < 2 * blockSize; ++y) {
					Pt a = Pt{
						.x = x,
						.y = y
					};
					bool bit = (Dist(doubledLine.a, doubledLine.b, a) >= 0.0);
					double value = bit ? +1.0 : -1.0;
					if (bit) {
						allClear = false;
					} else {
						allSet = false;
					}
					zoomIn[x + (2 * blockSize * y)] = value;
				}
			}
			double total = 0.0;
			for (int x = 0; x < blockSize; ++x) {
				for (int y = 0; y < blockSize; ++y) {
					double weight = 0.0;
					double tot = 0.0;
					for (int dx = -kernelHalfWidth; dx <= kernelHalfWidth; ++dx) {
						int u = std::clamp(2 * x + dx, 0, static_cast<int>(2 * blockSize) - 1);
						for (int dy = -kernelHalfWidth; dy <= kernelHalfWidth; ++dy) {
							int v = std::clamp(2 * y + dy, 0, static_cast<int>(2 * blockSize) - 1);
							double coeff = gaussian[(dx + kernelHalfWidth) + kernelWidth * (dy + kernelHalfWidth)];
							weight += coeff;
							tot += coeff * zoomIn[u + 2 * blockSize * v];
						}
					}
					double value = tot / weight;
					total += value;
					dictionary[i][x + blockSize * y] = value;
				}
			}
			double mean = total / static_cast<double>(blockSize * blockSize);
			double sumsq = 0.0;
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				double value = dictionary[i][j];
				if (!(allSet || allClear)) {
					value -= mean;
					dictionary[i][j] = value;
				}
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

	Eigen::MatrixXf createSegmentDictionaryFast(size_t blockSize, const std::vector<Line>& basisSet) {
		Eigen::MatrixXf dictionary(basisSet.size(), blockSize * blockSize);
		const double sigma = 1.0;
		const int kernelWidth = static_cast<int>(1 + sigma * 6);
		const int kernelHalfWidth = kernelWidth / 2;
		Eigen::MatrixXf gaussian(kernelWidth, kernelWidth);
		for (int dx = -kernelHalfWidth; dx <= +kernelHalfWidth; ++dx) {
			for (int dy = -kernelHalfWidth; dy <= +kernelHalfWidth; ++dy) {
				gaussian(dx + kernelHalfWidth, dy + kernelHalfWidth) = static_cast<float>(exp(-static_cast<double>(dx * dx + dy * dy) / (2.0 * sigma * sigma)));
			}
		}
		size_t i = 0;
		for (const Line separator : basisSet) {
			bool allSet = true;
			bool allClear = true;
			Eigen::VectorXf zoomIn(4 * blockSize * blockSize);
			Line doubledLine = Line{
				.a = Pt {
					.x = separator.a.x * 2,
					.y = separator.a.y * 2
				},
				.b = Pt {
					.x = separator.b.x * 2,
					.y = separator.b.y * 2
				}
			};
			for (int x = 0; x < 2 * blockSize; ++x) {
				for (int y = 0; y < 2 * blockSize; ++y) {
					Pt a = Pt{
						.x = x,
						.y = y
					};
					bool bit = (Dist(doubledLine.a, doubledLine.b, a) >= 0.0);
					float value = bit ? +1.0f : -1.0f;
					if (bit) {
						allClear = false;
					} else {
						allSet = false;
					}
					zoomIn[x + (2 * blockSize * y)] = value;
				}
			}
			float total = 0.0;
			for (int x = 0; x < blockSize; ++x) {
				for (int y = 0; y < blockSize; ++y) {
					float weight = 0.0;
					float tot = 0.0;
					for (int dx = -kernelHalfWidth; dx <= kernelHalfWidth; ++dx) {
						int u = std::clamp(2 * x + dx, 0, static_cast<int>(2 * blockSize) - 1);
						for (int dy = -kernelHalfWidth; dy <= kernelHalfWidth; ++dy) {
							int v = std::clamp(2 * y + dy, 0, static_cast<int>(2 * blockSize) - 1);
							float coeff = gaussian(dx + kernelHalfWidth, dy + kernelHalfWidth);
							weight += coeff;
							tot += coeff * zoomIn[u + 2 * blockSize * v];
						}
					}
					float value = tot / weight;
					total += value;
					dictionary(i, x + blockSize * y) = value;
				}
			}
			float mean = total / static_cast<float>(blockSize * blockSize);
			float sumsq = 0.0;
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				float value = dictionary(i, j);
				if (!(allSet || allClear)) {
					value -= mean;
					dictionary(i, j) = value;
				}
				sumsq += value * value;
			}
			float norm = sqrtf(sumsq);
			for (size_t j = 0; j < blockSize * blockSize; ++j) {
				float value = dictionary(i, j);
				if (sumsq != 0.0) {
					dictionary(i, j) = value / norm;
				}
			}
			++i;
		}
		return dictionary;
	}

	math::Matrix createIntraSegmentDictionary(size_t blockSize, const Line& segmentMask, CovarianceModel model) {
		std::vector<Pt> map1;
		std::vector<Pt> map2;
		for (int x = 0; x < blockSize; ++x) {
			for (int y = 0; y < blockSize; ++y) {
				Pt a = Pt{
					.x = x,
					.y = y
				};
				double dist = Dist(segmentMask.a, segmentMask.b, a);
				if (dist >= 0.0) {
					map1.push_back(a);
				} else {
					map2.push_back(a);
				}
			}
		}
		math::Matrix part1(map1.size(), map1.size());
		for (size_t i = 0; i < map1.size(); ++i) {
			const Pt& a = map1[i];
			for (size_t j = 0; j < map1.size(); ++j) {
				const Pt& b = map1[j];
				double dx = static_cast<double>(a.x - b.x);
				double dy = static_cast<double>(a.y - b.y);
				part1[i][j] = model(dx, dy);
			}
		}
		math::Matrix part2(map2.size(), map2.size());
		for (size_t i = 0; i < map2.size(); ++i) {
			const Pt& a = map2[i];
			for (size_t j = 0; j < map2.size(); ++j) {
				const Pt& b = map2[j];
				double dx = static_cast<double>(a.x - b.x);
				double dy = static_cast<double>(a.y - b.y);
				part2[i][j] = model(dx, dy);
			}
		}

		std::vector<math::Vector> basis1 = createBasis(part1);
		std::vector<math::Vector> basis2 = createBasis(part2);
		size_t part1Count = std::max(0, static_cast<int>(map1.size()) - 1);
		size_t part2Count = std::max(0, static_cast<int>(map2.size()) - 1);

		math::Matrix dictionary(part1Count + part2Count, blockSize * blockSize);
		math::Vector temp(blockSize * blockSize);
		size_t offset = 0;
		for (size_t i = 0; i < std::max(part1Count, part2Count); ++i) {
			if (i < part1Count) {
				memset(temp.Data(), 0, sizeof(double) * temp.Length());
				double meanTotal = 0.0;
				for (size_t j = 0; j < map1.size(); ++j) {
					const Pt& a = map1[j];
					double value = basis1[i + 1][j];
					meanTotal += value;
					temp[a.x + blockSize * a.y] = value;
				}
				meanTotal /= static_cast<double>(map1.size());
				double sumSq = 0.0;
				for (size_t j = 0; j < map1.size(); ++j) {
					const Pt& a = map1[j];
					double value = temp[a.x + blockSize * a.y];
					value -= meanTotal;
					sumSq += value * value;
					temp[a.x + blockSize * a.y] = value;
				}
				sumSq = sqrt(sumSq);
				for (size_t j = 0; j < blockSize * blockSize; ++j) {
					double value = temp[j];
					value /= sumSq;
					dictionary[offset][j] = value;
				}
				++offset;
			}
			if (i < part2Count) {
				double meanTotal = 0.0;
				memset(temp.Data(), 0, sizeof(double) * temp.Length());
				for (size_t j = 0; j < map2.size(); ++j) {
					const Pt& a = map2[j];
					double value = basis2[i + 1][j];
					meanTotal += value;
					temp[a.x + blockSize * a.y] = value;
				}
				meanTotal /= static_cast<double>(map2.size());
				double sumSq = 0.0;
				for (size_t j = 0; j < map2.size(); ++j) {
					const Pt& a = map2[j];
					double value = temp[a.x + blockSize * a.y];
					value -= meanTotal;
					sumSq += value * value;
					temp[a.x + blockSize * a.y] = value;
				}
				sumSq = sqrt(sumSq);
				for (size_t j = 0; j < blockSize * blockSize; ++j) {
					double value = temp[j];
					if (sumSq != 0.0) {
						value /= sumSq;
					}
					dictionary[offset][j] = value;
				}
				++offset;
			}
		}
		return dictionary;
	}

	Eigen::MatrixXf createIntraSegmentDictionaryFast(size_t blockSize, const Line& segmentMask, CovarianceModel model) {
		std::vector<Pt> map1;
		std::vector<Pt> map2;
		for (int x = 0; x < blockSize; ++x) {
			for (int y = 0; y < blockSize; ++y) {
				Pt a = Pt{
					.x = x,
					.y = y
				};
				double dist = Dist(segmentMask.a, segmentMask.b, a);
				if (dist >= 0.0) {
					map1.push_back(a);
				} else {
					map2.push_back(a);
				}
			}
		}
		Eigen::MatrixXf part1(map1.size(), map1.size());
		for (size_t i = 0; i < map1.size(); ++i) {
			const Pt& a = map1[i];
			for (size_t j = 0; j < map1.size(); ++j) {
				const Pt& b = map1[j];
				double dx = static_cast<double>(a.x - b.x);
				double dy = static_cast<double>(a.y - b.y);
				part1(i, j) = static_cast<float>(model(dx, dy));
			}
		}
		Eigen::MatrixXf part2(map2.size(), map2.size());
		for (size_t i = 0; i < map2.size(); ++i) {
			const Pt& a = map2[i];
			for (size_t j = 0; j < map2.size(); ++j) {
				const Pt& b = map2[j];
				double dx = static_cast<double>(a.x - b.x);
				double dy = static_cast<double>(a.y - b.y);
				part2(i, j) = static_cast<float>(model(dx, dy));
			}
		}

		std::vector<Eigen::VectorXf> basis1 = createBasisFast(part1);
		std::vector<Eigen::VectorXf> basis2 = createBasisFast(part2);
		size_t part1Count = std::max(0, static_cast<int>(map1.size()) - 1);
		size_t part2Count = std::max(0, static_cast<int>(map2.size()) - 1);

		Eigen::MatrixXf dictionary(part1Count + part2Count, blockSize * blockSize);
		Eigen::VectorXf temp(blockSize * blockSize);
		size_t offset = 0;
		for (size_t i = 0; i < std::max(part1Count, part2Count); ++i) {
			if (i < part1Count) {
				temp.setConstant(0.0);
				float meanTotal = 0.0;
				for (size_t j = 0; j < map1.size(); ++j) {
					const Pt& a = map1[j];
					float value = basis1[i + 1][j];
					meanTotal += value;
					temp[a.x + blockSize * a.y] = value;
				}
				meanTotal /= static_cast<float>(map1.size());
				float sumSq = 0.0;
				for (size_t j = 0; j < map1.size(); ++j) {
					const Pt& a = map1[j];
					float value = temp[a.x + blockSize * a.y];
					value -= meanTotal;
					sumSq += value * value;
					temp[a.x + blockSize * a.y] = value;
				}
				sumSq = sqrt(sumSq);
				for (size_t j = 0; j < blockSize * blockSize; ++j) {
					float value = temp[j];
					value /= sumSq;
					dictionary(offset, j) = value;
				}
				++offset;
			}
			if (i < part2Count) {
				float meanTotal = 0.0;
				temp.setConstant(0.0);
				for (size_t j = 0; j < map2.size(); ++j) {
					const Pt& a = map2[j];
					float value = basis2[i + 1](j);
					meanTotal += value;
					temp[a.x + blockSize * a.y] = value;
				}
				meanTotal /= static_cast<float>(map2.size());
				float sumSq = 0.0;
				for (size_t j = 0; j < map2.size(); ++j) {
					const Pt& a = map2[j];
					float value = temp[a.x + blockSize * a.y];
					value -= meanTotal;
					sumSq += value * value;
					temp[a.x + blockSize * a.y] = value;
				}
				sumSq = sqrt(sumSq);
				for (size_t j = 0; j < blockSize * blockSize; ++j) {
					float value = temp[j];
					if (sumSq != 0.0) {
						value /= sumSq;
					}
					dictionary(offset, j) = value;
				}
				++offset;
			}
		}
		return dictionary;
	}
}