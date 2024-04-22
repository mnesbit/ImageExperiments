#ifndef _BASIS_SET_H_
#define _BASIS_SET_H_

#include "framework.h"

#include <vector>
#include <functional>
#include "../../SimpleMatrix/inc/mathmatrix.h"
#include "../../SimpleMatrix/inc/mathvector.h"

namespace basis {
	typedef std::function<double(double, double)> CovarianceModel;

	double covariancePredictionModelY(double dx, double dy);

	double covariancePredictionModelU(double dx, double dy);

	double covariancePredictionModelV(double dx, double dy);

	std::vector<math::Vector> createBasis(math::Matrix& covariance);

	std::vector<std::vector<double> > distinctLineShapes(size_t blockSize);

	math::Matrix createSegmentDictionary(size_t blockSize, const std::vector<std::vector<double> >& basisSet);

	math::Matrix createIntraSegmentDictionary(size_t blockSize, const std::vector<double>& segmentMask, CovarianceModel model);
}


#endif //_BASIS_SET_H_