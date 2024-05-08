#include "../inc/MatchingPursuit.h"
#include "../inc/BitBuffer.h"
#include <algorithm>

namespace matching {

	static double Select(const math::Vector& residual, const math::Matrix& dictionary, math::Vector& selected, int& index)
	{
		index = -1;
		double bestCoeff = 0.0;
		math::Vector proj = Multiply(dictionary, residual);
		double* pData = proj.Data();
		size_t N = proj.Length();
		for (size_t i = 0; i < N; ++i) {
			double value = *pData++;
			if (abs(value) > abs(bestCoeff)) {
				bestCoeff = value;
				index = static_cast<int>(i);
			}
		}
		if (index != -1) {
			selected = dictionary.GetRow(index);
		}
		return bestCoeff;
	}

	static float SelectFast(const Eigen::VectorXf& residual, const Eigen::MatrixXf& dictionary,Eigen::VectorXf& selected, int& index)
	{
		const Eigen::VectorXf proj = dictionary * residual;
		index = -1;
		proj.cwiseAbs().maxCoeff(&index);
		if (index != -1) {
			selected = dictionary.row(index);
			return proj[index];
		}
		return 0.0;
	}

	int CalcMPDynamic(int K, const double quantization[], std::vector<BasisChoice>& results, const math::Vector& input, const DynamicDictionaryFunction& dynamicDictionary)
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
			coeff = quantization[reps] * static_cast<double>(quantized);
			if (quantized == 0)
			{
				return reps;
			}
			newEntry.Scale(coeff);
			residual.Subtract(newEntry);
		}
		return K;
	}

	int CalcMPDynamicFast(int K, const Eigen::VectorXf& quantization, std::vector<BasisChoice>& results, const Eigen::VectorXf& input, const DynamicDictionaryFunctionFast& dynamicDictionary)
	{
		size_t N = input.size();
		Eigen::VectorXf residual(input);
		Eigen::VectorXf newEntry(N);
		int prevBasis = 0;
		for (int reps = 0; reps < K; ++reps) {
			const Eigen::MatrixXf dictionary = dynamicDictionary(reps, results);
			int basisId;
			float coeff = SelectFast(residual, dictionary, newEntry, basisId);
			if (basisId < 0) {
				results[reps].deltaId = 0;
				results[reps].intCoeff = 0;
				return reps;
			}
			if (reps > 0) {
				int diff = basisId - prevBasis;
				results[reps].deltaId = static_cast<unsigned short>(bitbuffer::zigzagEncode(diff));
			} else {
				results[reps].deltaId = static_cast<unsigned short>(basisId);
			}
			prevBasis = basisId;
			int quantized = static_cast<int>(round(coeff / quantization[reps]));
			results[reps].intCoeff = static_cast<unsigned short>(bitbuffer::zigzagEncode(quantized));
			coeff = quantization[reps] * static_cast<float>(quantized);
			if (quantized == 0) {
				return reps;
			}
			residual -= (coeff * newEntry);
		}
		return K;
	}

	math::Vector FromCoeffsDynamic(int K, const double quantization[], const std::vector<BasisChoice>& coeffs, const DynamicDictionaryFunction& dynamicDictionary)
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

	Eigen::VectorXf FromCoeffsDynamicFast(int K, const Eigen::VectorXf& quantization, const std::vector<BasisChoice>& coeffs, const DynamicDictionaryFunctionFast& dynamicDictionary)
	{
		Eigen::MatrixXf dictionary = dynamicDictionary(K, coeffs);
		size_t N = dictionary.cols();
		Eigen::VectorXf results(N);
		results.setConstant(0.0);
		int choice = 0;
		for (int i = 0; i < K; ++i) {
			if (i > 0) {
				choice = choice + bitbuffer::zigzagDecode(static_cast<unsigned int>(coeffs[i].deltaId));
			} else {
				choice = static_cast<int>(coeffs[0].deltaId);
			}
			float coeff = quantization[i] * static_cast<float>(bitbuffer::zigzagDecode(static_cast<unsigned int>(coeffs[i].intCoeff)));
			results += coeff * dictionary.row(choice);
		}
		return results;
	}

}
