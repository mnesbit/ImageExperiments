#ifndef _MATCHING_PURSUIT_H_
#define _MATCHING_PURSUIT_H_

#include <functional>
#include "../SimpleMatrix/inc/mathmatrix.h"
#include "../SimpleMatrix/inc/mathvector.h"

namespace matching {

	typedef struct BasisChoice_t
	{
		unsigned short deltaId;
		unsigned short intCoeff;
	} BasisChoice;

	typedef std::function <math::Matrix(int, const std::vector<BasisChoice>)> DynamicDictionaryFunction;

	int CalcMPDynamic(int K, const double quantization[], std::vector<BasisChoice>& results, const math::Vector& input, DynamicDictionaryFunction dynamicDictionary);

	math::Vector FromCoeffsDynamic(int K, const double quantization[], const std::vector<BasisChoice>& coeffs, DynamicDictionaryFunction dynamicDictionary);
}

#endif //_MATCHING_PURSUIT_H_