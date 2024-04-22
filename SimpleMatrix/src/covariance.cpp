#include "../inc/covariance.h"

namespace math {

void Stat::update(double val) {
	if (N == 0) {
		N = 1.0;
		min = val;
		max = val;
		mean = val;
		sumSq = 0.0;
	} else {
		if (val < min) min = val;
		if (val > max) max = val;
		N += 1.0;
		double delta = val - mean;
		mean += delta / N;
		double delta2 = val - mean;
		sumSq += delta * delta2;
	}
}

double Stat::sampleVariance() const {
	return sumSq / (N - 1.0);
}

CovarianceCalculator::CovarianceCalculator()
	: m_N(0),
	m_means(),
	m_covariance()
{
}

CovarianceCalculator::CovarianceCalculator(const CovarianceCalculator& calculator)
	: m_N(calculator.m_N),
	m_means(calculator.m_means),
	m_covariance(calculator.m_covariance)
{

}

CovarianceCalculator::~CovarianceCalculator() {
	m_means.Reset();
	m_covariance.Reset();
}

CovarianceCalculator& CovarianceCalculator::operator=(const CovarianceCalculator& calculator) {
	if (this != &calculator)
	{
		m_N = calculator.m_N;
		m_means = calculator.m_means;
		m_covariance = calculator.m_covariance;
	}
	return *this;
}

void CovarianceCalculator::Update(Vector input) {
	if (m_N == 0) {
		m_means = Vector(input);
		m_covariance = Matrix(m_means.Length());
		m_N = 1;
		return;
	}
	if (m_means.Length() != input.Length()) {
		throw new std::range_error("mismatched sizes");
	}
	double n = static_cast<double>(++m_N);
	Vector diff(input);
	diff.Subtract(m_means);
	Vector scaledDiff(diff);
	scaledDiff.Scale(1.0 / n);
	m_means.Add(scaledDiff);
	for (int x = 0; x < m_means.Length(); ++x) {
		for (int y = 0; y < m_means.Length(); ++y) {
			m_covariance[x][y] = m_covariance[x][y] + (diff[x] * (input[y] - m_means[y]));
		}
	}
}


}