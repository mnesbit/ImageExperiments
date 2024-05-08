#ifndef _BASIS_SET_H_
#define _BASIS_SET_H_

#include "framework.h"

#include <vector>
#include <functional>
#include "../../SimpleMatrix/inc/mathmatrix.h"
#include "../../SimpleMatrix/inc/mathvector.h"
#include "../../eigen/Eigen/core"

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
	Eigen::MatrixXf createDCTDictionaryFast(int blockSize);

	math::Matrix createKLTDictionary(int blockSize, CovarianceModel model);
	Eigen::MatrixXf createKLTDictionaryFast(int blockSize, CovarianceModel model);

	std::vector<math::Vector> createBasis(const math::Matrix& covariance);
	std::vector<Eigen::VectorXf> createBasisFast(const Eigen::MatrixXf& covariance);

	math::Matrix createDictionary(const math::Matrix& covariance);
	Eigen::MatrixXf createDictionaryFast(const Eigen::MatrixXf& covariance);

	std::vector<Line> distinctLineShapes(size_t blockSize);

	math::Matrix createSegmentDictionary(size_t blockSize, const std::vector<Line>& basisSet);
	Eigen::MatrixXf createSegmentDictionaryFast(size_t blockSize, const std::vector<Line>& basisSet);

	math::Matrix createIntraSegmentDictionary(size_t blockSize, const Line& segmentMask, CovarianceModel model);
	Eigen::MatrixXf createIntraSegmentDictionaryFast(size_t blockSize, const Line& segmentMask, CovarianceModel model);
}


#endif //_BASIS_SET_H_