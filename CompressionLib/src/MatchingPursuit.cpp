#include "../inc/MatchingPursuit.h"
#include "../inc/BitBuffer.h"
#include <algorithm>

namespace matching {

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

	int CalcMPDynamic(int K, const double quantization[], std::vector<BasisChoice>& results, const math::Vector& input, DynamicDictionaryFunction dynamicDictionary)
	{
		size_t N = input.Length();
		math::Vector residual(input);
		math::Vector newEntry(N);
		int prevBasis = 0;
		for (int reps = 0; reps < K; ++reps)
		{
			math::Matrix dictionary = dynamicDictionary(reps, results);
			int basisId;
			double coeff = Select(residual, dictionary, newEntry, basisId);
			if (basisId < 0) {
				results[reps].deltaId = 0;
				results[reps].intCoeff = 0;
				return reps;
			}
			if (reps > 0) {
				int diff = basisId - prevBasis;
				results[reps].deltaId = static_cast<unsigned short>(bitbuffer::zigzagEncode(diff));
			}
			else {
				results[reps].deltaId = static_cast<unsigned short>(basisId);
			}
			prevBasis = basisId;
			int quantized = static_cast<int>(round(coeff / quantization[reps]));
			results[reps].intCoeff = static_cast<unsigned short>(bitbuffer::zigzagEncode(quantized));
			coeff = quantization[reps] * static_cast<double>(bitbuffer::zigzagDecode(results[reps].intCoeff));
			if (quantized == 0)
			{
				return reps;
			}
			newEntry.Scale(coeff);
			residual.Subtract(newEntry);
		}
		return K;
	}

	math::Vector FromCoeffsDynamic(int K, const double quantization[], const std::vector<BasisChoice>& coeffs, DynamicDictionaryFunction dynamicDictionary)
	{
		math::Matrix dictionary = dynamicDictionary(K, coeffs);
		size_t N = dictionary.Columns();
		math::Vector results(dictionary.Columns());
		int choice = 0;
		for (int i = 0; i < K; ++i) {
			if (i > 0) {
				choice = choice + bitbuffer::zigzagDecode(static_cast<unsigned int>(coeffs[i].deltaId));
			} else {
				choice = static_cast<int>(coeffs[0].deltaId);
			}
			double coeff = quantization[i] * static_cast<double>(bitbuffer::zigzagDecode(static_cast<unsigned int>(coeffs[i].intCoeff)));
			const math::MatrixRow basis = dictionary[choice];
			for (int j = 0; j < N; ++j) {
				results[j] += basis[j] * coeff;
			}
		}
		return results;
	}
}
