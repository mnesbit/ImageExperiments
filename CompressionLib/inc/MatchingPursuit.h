#ifndef _MATCHING_PURSUIT_H_
#define _MATCHING_PURSUIT_H_

#include "framework.h"

#include <functional>
#include "../../SimpleMatrix/inc/mathmatrix.h"
#include "../../SimpleMatrix/inc/mathvector.h"
#include "../../eigen/Eigen/core"

namespace matching {

	typedef struct BasisChoice_t
	{
		unsigned short deltaId;
		unsigned short intCoeff;
	} BasisChoice;

	typedef std::function <math::Matrix(int, const std::vector<BasisChoice>&)> DynamicDictionaryFunction;
	typedef std::function <Eigen::MatrixXf(int, const std::vector<BasisChoice>&)> DynamicDictionaryFunctionFast;

	int CalcMPDynamic(int K, const double quantization[], std::vector<BasisChoice>& results, const math::Vector& input, const DynamicDictionaryFunction& dynamicDictionary);
	int CalcMPDynamicFast(int K, const Eigen::VectorXf& quantization, std::vector<BasisChoice>& results, const Eigen::VectorXf& input, const DynamicDictionaryFunctionFast& dynamicDictionary);

	math::Vector FromCoeffsDynamic(int K, const double quantization[], const std::vector<BasisChoice>& coeffs, const DynamicDictionaryFunction& dynamicDictionary);
	Eigen::VectorXf FromCoeffsDynamicFast(int K, const Eigen::VectorXf& quantization, const std::vector<BasisChoice>& coeffs, const DynamicDictionaryFunctionFast& dynamicDictionary);
}

#endif //_MATCHING_PURSUIT_H_