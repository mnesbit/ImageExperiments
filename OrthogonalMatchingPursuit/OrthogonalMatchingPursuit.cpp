#include <iostream>
#include <vector>
#include <set>
#include <format>
#include <algorithm>
#include <functional>
#include "../ImageHelper/inc/imgloader.h"
#include "../SimpleMatrix/inc/covariance.h"
#include "../SimpleMatrix/inc/symmeigen.h"
#include "../ImageHelper/inc/misc.h"

using namespace img;

typedef struct BasisChoice_Tag
{
	int basisID;
	double coeff;
}
BasisChoice;


typedef std::function <math::Matrix(int, const BasisChoice[])> DynamicDictionaryFunction;

double Select(const math::Vector& residual, const math::Matrix& dictionary, math::Vector& selected, int& index)
{
	index = -1;
	double bestCoeff = 0.0;
	math::Vector proj = Multiply(dictionary, residual);
	double* pData = proj.Data();
	size_t N = proj.Length();
	for (size_t i = 0; i < N; ++i)
	{
		double value = *pData++;
		if (abs(value) > abs(bestCoeff))
		{
			bestCoeff = value;
			index = static_cast<int>(i);
		}
	}
	if (index != -1)
	{
		selected = dictionary.GetRow(index);
	}
	return bestCoeff;
}

double CalcOMPStatic(int K, BasisChoice results[], const math::Vector& input, const math::Matrix& dictionary)
{
	size_t N = input.Length();
	math::Vector residual(input);
	math::Vector newEntry(N);
	math::Matrix Q(N, K);
	math::Matrix R(K, K);
	math::Vector z(K);
	math::Vector X(residual);
	for (int reps = 0; reps < K; ++reps)
	{
		results[reps].coeff = Select(residual, dictionary, newEntry, results[reps].basisID);
		if (abs(results[reps].coeff) < 1E-10)
		{
			for (int i = reps; i < K; ++i)
			{
				results[i].basisID = -1;
				results[i].coeff = 0.0;
			}
			K = reps;
			if (reps == 0)return 0.0;
			break;
		}
		for (int j = 0; j < reps; ++j)
		{
			math::Vector prevBasis = Q.GetColumn(j);
			double nonOrth = newEntry.DotProduct(prevBasis);
			R[j][reps] = nonOrth;
			prevBasis.Scale(nonOrth);
			newEntry.Subtract(prevBasis);
		}
		double mag = newEntry.Magnitude();
		R[reps][reps] = mag;
		newEntry.Scale(1.0 / mag);
		for (int j = 0; j < N; ++j)
		{
			Q[j][reps] = newEntry[j];
		}
		double coeff = newEntry.DotProduct(X);
		z[reps] = coeff;
		newEntry.Scale(coeff);
		residual.Subtract(newEntry);
	}
	z[K - 1] = z[K - 1] / R[K - 1][K - 1];
	results[K - 1].coeff = z[K - 1];
	for (int i = K - 2; i >= 0; --i)
	{
		double tot = 0.0;
		for (int j = i + 1; j < K; ++j)
		{
			tot += R[i][j] * z[j];
		}
		z[i] = (z[i] - tot) / R[i][i];
		results[i].coeff = z[i];
	}
	return residual.Magnitude();
}

double CalcOMPDynamic(int K, BasisChoice results[], const math::Vector& input, DynamicDictionaryFunction dynamicDictionary)
{
	size_t N = input.Length();
	math::Vector residual(input);
	math::Vector newEntry(N);
	math::Matrix Q(N, K);
	math::Matrix R(K, K);
	math::Vector z(K);
	math::Vector X(residual);
	for (int reps = 0; reps < K; ++reps)
	{
		math::Matrix dictionary = dynamicDictionary(reps, results);
		results[reps].coeff = Select(residual, dictionary, newEntry, results[reps].basisID);
		if (abs(results[reps].coeff) < 1E-10)
		{
			for (int i = reps; i < K; ++i)
			{
				results[i].basisID = -1;
				results[i].coeff = 0.0;
			}
			K = reps;
			if (reps == 0)return 0.0;
			break;
		}
		for (int j = 0; j < reps; ++j)
		{
			math::Vector prevBasis = Q.GetColumn(j);
			double nonOrth = newEntry.DotProduct(prevBasis);
			R[j][reps] = nonOrth;
			prevBasis.Scale(nonOrth);
			newEntry.Subtract(prevBasis);
		}
		double mag = newEntry.Magnitude();
		R[reps][reps] = mag;
		newEntry.Scale(1.0 / mag);
		for (int j = 0; j < N; ++j)
		{
			Q[j][reps] = newEntry[j];
		}
		double coeff = newEntry.DotProduct(X);
		z[reps] = coeff;
		newEntry.Scale(coeff);
		residual.Subtract(newEntry);
	}
	z[K - 1] = z[K - 1] / R[K - 1][K - 1];
	results[K - 1].coeff = z[K - 1];
	for (int i = K - 2; i >= 0; --i)
	{
		double tot = 0.0;
		for (int j = i + 1; j < K; ++j)
		{
			tot += R[i][j] * z[j];
		}
		z[i] = (z[i] - tot) / R[i][i];
		results[i].coeff = z[i];
	}
	return residual.Magnitude();
}

math::Vector FromCoeffsStatic(int K, const BasisChoice coeffs[], const math::Matrix& dictionary,const math::Vector original, double squaredError[], double summedMagnitude[]) {
	size_t N = dictionary.Columns();
	math::Vector results(dictionary.Columns());
	for (int i = 0; i < K; ++i) {
		double total = 0.0;
		if (coeffs[i].basisID >= 0 && coeffs[i].coeff != 0.0) {
			const math::MatrixRow basis = dictionary[coeffs[i].basisID];
			for (int j = 0; j < N; ++j) {
				int mag;
				frexp(coeffs[i].coeff, &mag);
				summedMagnitude[i] += static_cast<double>(mag);
				results[j] += basis[j] * static_cast<float>(coeffs[i].coeff); //crop coefficient to float to simulate some quantization
				double error = results[j] - original[j];
				total += error * error;
			}
		} else {
			break;
		}
		squaredError[i] += total;
	}
	return results;
}

math::Vector FromCoeffsDynamic(int K, const BasisChoice coeffs[], DynamicDictionaryFunction dynamicDictionary, const math::Vector original, double squaredError[], double summedMagnitude[]) {
	math::Matrix dictionary = dynamicDictionary(K, coeffs);
	size_t N = dictionary.Columns();
	math::Vector results(dictionary.Columns());
	for (int i = 0; i < K; ++i) {
		double total = 0.0;
		if (coeffs[i].basisID >= 0 && coeffs[i].coeff != 0.0) {
			const math::MatrixRow basis = dictionary[coeffs[i].basisID];
			for (int j = 0; j < N; ++j) {
				int mag;
				frexp(coeffs[i].coeff, &mag);
				summedMagnitude[i] += static_cast<double>(mag);
				results[j] += basis[j] * static_cast<float>(coeffs[i].coeff); //crop coefficient to float to simulate some quantization
				double error = results[j] - original[j];
				total += error * error;
			}
		} else {
			break;
		}
		squaredError[i] += total;
	}
	return results;
}

math::Matrix createDCTDictionary(int blockSize) {
	math::Matrix dctDict(blockSize * blockSize, blockSize * blockSize);
	for (int u0 = 0; u0 < blockSize; ++u0)
	{
		for (int v0 = 0; v0 < blockSize; ++v0)
		{
			int index = (u0 * blockSize) + v0;
			for (int x = 0; x < blockSize; ++x)
			{
				for (int y = 0; y < blockSize; ++y)
				{
					double sum = cos(((2 * x + 1) / (2.0 * blockSize)) * u0 * M_PI) * cos(((2 * y + 1) / (2.0 * blockSize)) * v0 * M_PI);
					if (u0 && v0)
					{
						sum = sum * sqrt(4.0 / (blockSize * blockSize));
					}
					else if (u0 || v0)
					{
						sum = sum * sqrt(2.0 / (blockSize * blockSize));
					}
					else
					{
						sum = 0.5 * sum * sqrt(4.0 / (blockSize * blockSize));
					}
					dctDict[index][x + (blockSize * y)] = sum;
				}
			}
		}
	}
	return dctDict;
}

// fully substituted version of optimised summedLaplaceModel from FitCovarianceModel experiments
double covariancePredictionModel(double dx, double dy) {
	return 3817.7299999999996 * exp(-1.48854e-05 * dx * dx + -1.7273e-05 * dy * dy)
		+ 657.8100000000001 * exp(-0.0436057 * abs(dx) + -0.050844400000000005 * abs(dy));
}

std::vector<math::Vector> createBasis(math::Matrix& part)
{
	math::SymmetricEigenDecomposition eigenDeco(part);
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

math::Matrix createKLTDictionary(int blockSize) {
	math::Matrix covarianceMatrix(blockSize * blockSize, blockSize * blockSize);
	for (int i = 0; i < blockSize * blockSize; ++i)
	{
		int xpos1 = i / blockSize;
		int ypos1 = i % blockSize;
		for (int j = 0; j < blockSize * blockSize; ++j)
		{
			int xpos2 = j / blockSize;
			int ypos2 = j % blockSize;
			covarianceMatrix[i][j] = covariancePredictionModel(xpos1 - xpos2, ypos1-ypos2);
		}
	}
	std::vector<math::Vector> sortedBasis = createBasis(covarianceMatrix);
	math::Matrix dictionary(sortedBasis.size(), blockSize * blockSize);
	for (int i = 0; i < sortedBasis.size(); ++i)
	{
		if (i == 0) { // ensure we have a flat DC basis
			for (int j = 0; j < blockSize * blockSize; ++j)
			{
				dictionary[i][j] = 1.0 / blockSize;
			}
		} else {
			for (int j = 0; j < blockSize * blockSize; ++j)
			{
				dictionary[i][j] = sortedBasis[i][j];
			}
		}
	}
	return dictionary;
}

typedef struct Pt_Tag
{
	int x, y;
} Pt;

static double Dist(Pt& line1, Pt& line2, Pt& point) {
	return static_cast<double>((line2.x - line1.x) * (line1.y - point.y) - (line1.x - point.x) * (line2.y - line1.y)) / sqrt((line2.x - line1.x) * (line2.x - line1.x) + (line2.y - line1.y) * (line2.y - line1.y));
}

struct CustomComparator {
	bool operator()(std::vector<double> a, std::vector<double> b) const
	{
		for (int i = 0; i < a.size(); ++i) {
			if (signbit(a[i]) != signbit(b[i])) {
				return signbit(a[i]);
			}
		}
		return false;
	}
};

std::vector<std::vector<double> > distinctLineShapes(int blockSize) {
	std::set<std::vector<double>, CustomComparator> basisSet;
	basisSet.insert(std::vector<double>(blockSize * blockSize, true)); //dc basis
	for (int side1 = -1; side1 < blockSize + 1; ++side1)
	{
		Pt s1 = { .x = side1, .y = 0 };
		for (int side2 = -1; side2 < blockSize + 1; ++side2)
		{
			Pt s2a = { .x = 0, .y = static_cast<int>(side2) };
			Pt s2b = { .x = static_cast<int>(side2), .y = blockSize - 1 };
			Pt s2c = { .x = blockSize - 1, .y = static_cast<int>(side2) };
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
	for (int side1 = -1; side1 < blockSize + 1; ++side1)
	{
		Pt s1 = { .x = blockSize - 1 , .y = side1 };
		for (int side2 = -1; side2 < blockSize + 1; ++side2)
		{
			Pt s2a = { .x = 0, .y = static_cast<int>(side2) };
			Pt s2b = { .x = static_cast<int>(side2), .y = blockSize - 1 };
			Pt s2c = { .x = 0, .y = side1 };
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

math::Matrix createSegmentDictionary(int blockSize, std::vector<std::vector<double> > basisSet) {
	math::Matrix dictionary(basisSet.size(), blockSize * blockSize);
	int i = 0;
	for (const std::vector<double>& bits: basisSet) {
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
		for (int j = 0; j < blockSize * blockSize; ++j) {
			double value = dictionary[i][j];
			if (!(allSet || allClear)) {
				value -= mean;
			}
			dictionary[i][j] = value;
			sumsq += value * value;
		}
		double norm = sqrt(sumsq);
		for (int j = 0; j < blockSize * blockSize; ++j) {
			double value = dictionary[i][j];
			if (sumsq != 0.0) {
				dictionary[i][j] = value / norm;
			}
		}
		++i;
	}
	return dictionary;
}

math::Matrix createSegmentDictionary(int blockSize) {
	std::vector<std::vector<double> > basisSet = distinctLineShapes(blockSize);
	return createSegmentDictionary(blockSize, basisSet);
}

void writeDictionaryToPNG(math::Matrix dictionary, std::string fileName) {
	size_t blockSize = static_cast<size_t>(sqrt(dictionary.Columns()));
	if (dictionary.Columns() != blockSize * blockSize) {
		throw std::invalid_argument("Dictionary should have square Column count");
	}
	size_t approxWidth = static_cast<size_t>(sqrt(dictionary.Rows()));
	size_t approxHeight = (dictionary.Rows() / approxWidth) + ((dictionary.Rows() % approxWidth == 0) ? 0 : 1);
	image<double>* basisPic = new image<double>(approxWidth * blockSize, approxHeight * blockSize);
	for (size_t i = 0; i < dictionary.Rows(); ++i) {
		size_t blockx = i % approxWidth;
		size_t blocky = i / approxWidth;
		for (size_t x = 0; x < blockSize; ++x) {
			for (size_t y = 0; y < blockSize; ++y) {
				imRef(basisPic, x + blockx * blockSize, y + blocky * blockSize) = dictionary[i][x + (y * blockSize)];
			}
		}
	}
	SaveImageGeneric(basisPic, fileName.c_str(), imgFormat::PNG);
	delete basisPic;
}

enum class Mode {
	DCT,
	KLT,
	SEG,
	SEG_KLT
};

int main(int argc, char* argv[])
{
	if (argc != 6 && argc != 7) {
		std::cout << "usage: OrthogonalMatchingPursuit.exe <BlockSize> <how many coefficients> <DCT|KLT|SEG|SEG_KLT> <input file> <output file> [<optional file to write basis as png image>]" << std::endl;
		return -1;
	}
	const int BlockSize = atoi(argv[1]);
	const int K = atoi(argv[2]);
	const std::string modeStr = argv[3];
	Mode mode;
	if (modeStr == "DCT") {
		mode = Mode::DCT;
		std::cout << std::format("using standard DCT basis calculated over block size of {} x {}", BlockSize, BlockSize) << std::endl;
	} else if (modeStr == "KLT") {
		mode = Mode::KLT;
		std::cout << std::format("using square Karhunen-Loeve basis derived from a covariance model and calculated over block size of {} x {}", BlockSize, BlockSize) << std::endl;
	} else if (modeStr == "SEG") {
		mode = Mode::SEG;
		std::cout << std::format("using a normalized basis derived from all distinct line cuts across the pixelated block of size {} x {}", BlockSize, BlockSize) << std::endl;
	} else if (modeStr == "SEG_KLT") {
		mode = Mode::SEG_KLT;
		std::cout << std::format("using a dynamic basis derived from all distinct line cuts across the pixelated block of size {} x {} and then introducing higher order KLT vectors over the two halves if the basic split is chosen", BlockSize, BlockSize) << std::endl;
	} else {
		std::cout << "mode not recognized: " << modeStr << " valid modes are DCT, KLT, SEG, or SEG_KLT case-sensitive" << std::endl;
		return -1;
	}
	std::cout << "Calculating basis vectors" << std::endl;
	math::Matrix baseDict;
	math::Matrix* detailBasis = nullptr;
	switch (mode) {
		case Mode::DCT:
			baseDict = createDCTDictionary(BlockSize);
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				writeDictionaryToPNG(baseDict, basisFileName);
			}
			break;
		case Mode::KLT:
			baseDict = createKLTDictionary(BlockSize);
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				writeDictionaryToPNG(baseDict, basisFileName);
			}
			break;
		case Mode::SEG:
			baseDict = createSegmentDictionary(BlockSize);
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				writeDictionaryToPNG(baseDict, basisFileName);
			}
			break;
		case Mode::SEG_KLT: {
			std::vector<std::vector<double> > basisSet = distinctLineShapes(BlockSize);
			baseDict = createSegmentDictionary(BlockSize, basisSet);
			size_t fullRowCount = baseDict.Rows();
			detailBasis = new math::Matrix[baseDict.Rows()];
			size_t id = 0;
			for (const std::vector<double>& basis : basisSet) {
				int count1 = 0;
				int count2 = 0;
				math::Matrix part1(BlockSize * BlockSize, BlockSize * BlockSize);
				math::Matrix part2(BlockSize * BlockSize, BlockSize * BlockSize);
				for (int x = 0; x < BlockSize * BlockSize; ++x) {
					int x1 = x % BlockSize;
					int y1 = x / BlockSize;
					if (basis[x] >= 0.0) ++count1;
					if (basis[x] <= 0.0) ++count2;
					for (int y = 0; y < BlockSize * BlockSize; ++y) {
						int x2 = y % BlockSize;
						int y2 = y / BlockSize;
						part1[x][y] = std::max(0.0, 0.5 + basis[x]) * std::max(0.0, 0.5 + basis[y]) * covariancePredictionModel(x2 - x1, y2 - y1);
						part2[x][y] = std::max(0.0, 0.5 - basis[x]) * std::max(0.0, 0.5 - basis[y]) * covariancePredictionModel(x2 - x1, y2 - y1);
					}
				}

				std::vector<math::Vector> basis1 = createBasis(part1);
				std::vector<math::Vector> basis2 = createBasis(part2);
				int part1Count = std::max(0, count1 - 2);
				int part2Count = std::max(0, count2 - 2);

				detailBasis[id] = math::Matrix(part1Count + part2Count, BlockSize * BlockSize);
				for (int i = 0; i < part1Count; ++i) {
					for (int j = 0; j < BlockSize * BlockSize; ++j) {
						detailBasis[id][i][j] = basis1[i + 1][j];
					}
				}
				for (int i = 0; i < part2Count; ++i) {
					for (int j = 0; j < BlockSize * BlockSize; ++j) {
						detailBasis[id][part1Count + i][j] = basis2[i + 1][j];
					}
				}
				fullRowCount += detailBasis[id].Rows();
				++id;
			}
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				math::Matrix fullBasis(fullRowCount, BlockSize * BlockSize);
				int offset = 0;
				for (int i = 0; i < baseDict.Rows(); ++i) {
					for (int j = 0; j < BlockSize * BlockSize; ++j) {
						fullBasis[offset][j] = baseDict[i][j];
					}
					++offset;
				}
				for (int subDict = 0; subDict < id; ++subDict) {
					for (int i = 0; i < detailBasis[subDict].Rows(); ++i) {
						for (int j = 0; j < BlockSize * BlockSize; ++j) {
							fullBasis[offset][j] = detailBasis[subDict][i][j];
						}
						++offset;
					}
				}
				writeDictionaryToPNG(fullBasis, basisFileName);
			}
		}
		break;
	}

	auto dynamic = [&] (int prevCoeffs, const BasisChoice prevChoices[])  -> math::Matrix {
		size_t atomCount = baseDict.Rows();
		std::vector<int> chosen;
		for (int i = 0; i < prevCoeffs; ++i) {
			int choice = prevChoices[i].basisID;
			if (choice >= 0
				&& choice < baseDict.Rows()
				&& std::find(chosen.cbegin(),chosen.cend(),choice) == chosen.cend()) {
				chosen.push_back(choice);
				atomCount += detailBasis[choice].Rows();
			}
		}
		math::Matrix dictionary(atomCount, BlockSize * BlockSize);
		int offset = 0;
		for (int i = 0; i < baseDict.Rows(); ++i) {
			for (int j = 0; j < BlockSize * BlockSize; ++j) {
				dictionary[offset][j] = baseDict[i][j];
			}
			++offset;
		}
		for (int choice : chosen) {
			math::Matrix& detail = detailBasis[choice];
			for (int i = 0; i < detail.Rows(); ++i) {
				for (int j = 0; j < BlockSize * BlockSize; ++j) {
					dictionary[offset][j] = detail[i][j];
				}
				++offset;
			}
		}
		return dictionary;
	};

	image<double>* imgIn = LoadImageGenericMono(argv[4]);
	std::cout << std::format("start processing image {} using largest {} coefficients.", argv[4], K) << std::endl;
	const size_t width = imgIn->width();
	const size_t height = imgIn->height();
	image<double>* imgOut = new image<double>(width, height, false);
	double* squaredError = new double[K];
	double* summedMagnitude = new double[K];
	for (int i = 0; i < K; ++i) {
		squaredError[i] = 0.0;
		summedMagnitude[i] = 0.0;
	}
	for (size_t x = 0; x < width; x += BlockSize) {
		std::cout << std::format("{} %", (100 * x) / width) << std::endl;
		for (size_t y = 0; y < height; y += BlockSize) {
			math::Vector block(BlockSize * BlockSize);
			for (size_t dx = 0; dx < BlockSize; ++dx) {
				size_t u = x + dx;
				for (size_t dy = 0; dy < BlockSize; ++dy) {
					size_t v = y + dy;
					if ((u < width) && (v < height))
					{
						block[dx + (BlockSize * dy)] = imRef(imgIn, u, v);
					}
					else
					{
						block[dx + (BlockSize * dy)] = 0.0;
					}
				}
			}
			BasisChoice* choices = new BasisChoice[K];
			math::Vector decoded;
			switch (mode) {
				case Mode::DCT:
				case Mode::KLT:
				case Mode::SEG: {
					CalcOMPStatic(K, choices, block, baseDict);
					decoded = FromCoeffsStatic(K, choices, baseDict, block, squaredError, summedMagnitude);
					break;
				}
				case Mode::SEG_KLT: {
					CalcOMPDynamic(K, choices, block, dynamic);
					decoded = FromCoeffsDynamic(K, choices, dynamic, block, squaredError, summedMagnitude);
					break;
				}
			}
			for (size_t dx = 0; dx < BlockSize; ++dx) {
				size_t u = x + dx;
				for (size_t dy = 0; dy < BlockSize; ++dy) {
					size_t v = y + dy;
					if ((u < width) && (v < height))
					{
						imRef(imgOut, u, v) = decoded[dx + (BlockSize * dy)];
					}
				}
			}
			delete[] choices;

		}
	}
	std::cout << std::format("Results of {} mode up to {} coefficients on {} x {} blocks", modeStr, K, BlockSize, BlockSize) << std::endl;
	for (int i = 0; i < K; ++i) {
		double mse = squaredError[i] / (width * height);
		double psnr = 20.0 * log10(255.0) - 10.0 * log10(mse);
		double averageMag = summedMagnitude[i] / (width * height);
		std::cout << std::format("PSNR at {} coefficients {} average magnitude {}", i + 1, psnr, averageMag) << std::endl;
	}
	SaveImageGeneric(imgOut, argv[5], imgFormat::PNG);
	delete imgIn;
	delete imgOut;
	delete[] squaredError;
	delete[] summedMagnitude;
	delete[] detailBasis;
    return 0;
}

