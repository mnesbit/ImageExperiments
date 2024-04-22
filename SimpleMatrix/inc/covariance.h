#ifndef _COVARIANCE_H_
#define _COVARIANCE_H_

#include "framework.h"

#include "mathvector.h"
#include "mathmatrix.h"

namespace math {
	struct Stat {
		double N;
		double min;
		double max;
		double mean;
		double sumSq;

		void update(double val);

		double sampleVariance() const;
	};


	class CovarianceCalculator {
	public:
		CovarianceCalculator();
		CovarianceCalculator(const CovarianceCalculator& calculator);
		~CovarianceCalculator();

		CovarianceCalculator& operator=(const CovarianceCalculator& calculator);

		void Update(Vector input);

		int Length() const { return m_N;  }

		Vector Means() const { return Vector(m_means);  }

		Matrix Covariance() const {
			Matrix cov(m_covariance);
			cov.Scale(1.0 / m_N);
			return cov;
		}
	private:
		int m_N;
		Vector m_means;
		Matrix m_covariance;
	};
}

#endif // _COVARIANCE_H_