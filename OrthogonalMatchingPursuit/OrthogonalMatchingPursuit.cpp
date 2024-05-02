#include <iostream>
#include <vector>
#include <set>
#include <format>
#include <algorithm>
#include <functional>
#include "../ImageHelper/inc/imgloader.h"
#include "../SimpleMatrix/inc/covariance.h"
#include "../SimpleMatrix/inc/symmeigen.h"
#include "../CompressionLib/inc/BasisSet.h"

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
	std::vector<math::Matrix> detailBasis;
	switch (mode) {
		case Mode::DCT:
			baseDict = basis::createDCTDictionary(BlockSize);
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				writeDictionaryToPNG(baseDict, basisFileName);
			}
			break;
		case Mode::KLT:
			baseDict = basis::createKLTDictionary(BlockSize,basis::covariancePredictionModelY);
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				writeDictionaryToPNG(baseDict, basisFileName);
			}
			break;
		case Mode::SEG:
		{
			std::vector<basis::Line> shapes = basis::distinctLineShapes(BlockSize);
			baseDict = basis::createSegmentDictionary(BlockSize, shapes);
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				writeDictionaryToPNG(baseDict, basisFileName);
			}
		}
			break;
		case Mode::SEG_KLT: {
			std::vector<basis::Line> shapes = basis::distinctLineShapes(BlockSize);
			baseDict = basis::createSegmentDictionary(BlockSize, shapes);
			size_t fullRowCount = baseDict.Rows();
			for (const basis::Line& shape: shapes) {
				math::Matrix detail = basis::createIntraSegmentDictionary(BlockSize, shape, basis::covariancePredictionModelY);
				fullRowCount += detail.Rows();
				detailBasis.push_back(std::move(detail));
			}
			if (argc == 7) {
				std::string basisFileName(argv[6]);
				math::Matrix fullBasis(fullRowCount, BlockSize * BlockSize);
				size_t offset = 0;
				for (size_t i = 0; i < baseDict.Rows(); ++i) {
					for (size_t j = 0; j < BlockSize * BlockSize; ++j) {
						fullBasis[offset][j] = baseDict[i][j];
					}
					++offset;
				}
				for (size_t subDict = 0; subDict < shapes.size(); ++subDict) {
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
    return 0;
}

