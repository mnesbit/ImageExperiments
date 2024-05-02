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

	typedef struct Pt_Tag
	{
		int x, y;
	} Pt;

	typedef struct Line_Tag
	{
		Pt a;
		Pt b;
	} Line;

	math::Matrix createDCTDictionary(int blockSize);

	math::Matrix createKLTDictionary(int blockSize, CovarianceModel model);

	std::vector<math::Vector> createBasis(math::Matrix& covariance);

	math::Matrix createDictionary(math::Matrix& covariance);

	std::vector<Line> distinctLineShapes(size_t blockSize);

	math::Matrix createSegmentDictionary(size_t blockSize, const std::vector<Line>& basisSet);

	math::Matrix createIntraSegmentDictionary(size_t blockSize, const Line& segmentMask, CovarianceModel model);
}


#endif //_BASIS_SET_H_