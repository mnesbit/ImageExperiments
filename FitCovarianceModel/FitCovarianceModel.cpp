#include <iostream>
#include <filesystem>
#include <random>
#include <limits>
#include <functional>
#include "../ImageHelper/inc/imgloader.h"
#include "../SimpleMatrix/inc/covariance.h"
#include "../ImageHelper/inc/misc.h"

using namespace img;

//final error: 9.17026e+10
//0.826782, 4.67431e-05, 0.146512, 0.725486, 0.0230279,
double linearModel(double params[5], double dx, double dy) {
	return 5000.0 * params[0] - params[1] * dx * dx - params[2] * dy * dy - 10.0 * params[3] * abs(dx) - 10.0 * params[4] * abs(dy);
}

//final error: 1.14187e+11
//0.813564, 0.368587,
double circularLaplaceModel(double params[2], double dx, double dy) {
	return 5000.0 * params[0] * exp(-0.0001 * params[1] * (dx * dx + dy * dy));
}

//final error: 1.12314e+11
//0.813581, 0.345163, 0.392403,
double ellipticLaplaceModel(double params[3], double dx, double dy) {
	return 5000.0 * params[0] * exp(-(0.0001 * params[1] * dx * dx + 0.0001 * params[2] * dy * dy));
}

//final error: 4.58244e+10
//0.846162, 0.210703,
double manhattanLaplaceModel(double params[2], double dx, double dy) {
	return 5000.0 * params[0] * exp(-0.01 * params[1] * (abs(dx) + abs(dy)));
}

//final error: 4.34405e+10
//0.846215, 0.197056, 0.224601,
double rectangularLaplaceModel(double params[3], double dx, double dy) {
	return 5000.0 * params[0] * exp(-(0.01 * params[1] * abs(dx) + 0.01 * params[2] * abs(dy)));
}

//final error: 4.34422e+10
//0.846182, 0, 5.03342e-06, 0.196322, 0.22491,
double combinedLaplaceModel(double params[5], double dx, double dy) {
	return 5000.0 * params[0] * exp(-(0.0001 * params[1] * dx * dx + 0.0001 * params[2] * dy * dy + 0.01 * params[3] * abs(dx) + 0.01 * params[4] * abs(dy)));
}

//final error: 5.78702e+09
//0.74988, 0.116278, 0.151871, 0.133772, 0.349443, 0.365449,
//after 1000 iterations : 3.99418e+09
//0.763546, 0.148854, 0.17273, 0.131562, 0.436057, 0.508444,
double summedLaplaceModel(double params[6], double dx, double dy) {
	return 5000.0 * params[0] * exp(-(0.0001 * params[1] * dx * dx + 0.0001 * params[2] * dy * dy))
		+ 5000.0 * params[3] * exp(-(0.1 * params[4] * abs(dx) + 0.1 * params[5] * abs(dy)));
}

// fully substituted version of optimised summedLaplaceModel
double finalModel(double dx, double dy) {
	return 3817.7299999999996 * exp(-1.48854e-05 * dx * dx + -1.7273e-05 * dy * dy)
		+ 657.8100000000001 * exp(-0.0436057 * abs(dx) + -0.050844400000000005 * abs(dy));
}

double modelError(math::Matrix& covar, int N, double params[]) {
	double error = 0.0;
	for (int x = 0; x < covar.Rows(); ++x) {
		int x1 = x % N;
		int y1 = x / N;
		math::MatrixRow row{ covar[x] };
		for (int y = 0; y < covar.Columns(); ++y) {
			int x2 = y % N;
			int y2 = y / N;
			double dx = static_cast<double>(x2 - x1);
			double dy = static_cast<double>(y2 - y1);
			double estimate = summedLaplaceModel(params, dx, dy);
			double actual = row[y];
			double diff = (actual - estimate);
			error += diff * diff;
		}
	}

	return error;
}

void particleSwarm(int dims, int particles, int iterations, std::function<double(double*, int)> modelCalc) {
	std::mt19937 rand;
	rand.seed(777);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	const double c1 = 2.0;
	const double c2 = 2.0;
	const double inertia = 0.5;
	double* positions = new double[particles * dims];
	double* bestPositions = new double[particles * dims];
	double* velocities = new double[particles * dims];
	double* bestValues = new double[particles];
	double* globalBestPoint = new double[dims];
	double* probe = new double[dims];
	double globalBest = 0.0;
	for (int i = 0; i < particles; ++i) {
		for (int j = 0; j < dims; ++j) {
			int offset = i * dims + j;
			positions[offset] = unif(rand);
			bestPositions[offset] = positions[offset];
			velocities[offset] = 1.0 - 2.0 * unif(rand);
		}
		bestValues[i] = modelCalc(positions + i * dims, dims);
		if ((i == 0) || (bestValues[i] < globalBest)) {
			globalBest = bestValues[i];
			for (int j = 0; j < dims; ++j) {
				globalBestPoint[j] = positions[i * dims + j];
			}
		}
	}
	std::cout << "initial error: " << globalBest << std::endl;
	for (int i = 0; i < dims; ++i) {
		std::cout << globalBestPoint[i] << ", ";
	}
	std::cout << std::endl;
	for (int iter = 0; iter < iterations; ++iter) {
		for (int i = 0; i < particles; i++)
		{
			if (unif(rand) < 0.05) //spawn a fresh agent to ensure innovation
			{
				bestValues[i] = std::numeric_limits<double>::max();
				for (int j = 0; j < dims; j++)
				{
					int offset = i * dims + j;
					positions[offset] = unif(rand);
					bestPositions[offset] = positions[offset];
					velocities[offset] = 1.0 - 2.0 * unif(rand);
				}
			}
			double rP = unif(rand);
			double rG = unif(rand);
			for (int j = 0; j < dims; j++)
			{
				int offset = i * dims + j;
				velocities[offset] = (inertia * velocities[offset]) + (c1 * rP * (bestPositions[offset] - positions[offset]));
				velocities[offset] += c2 * rG * (globalBestPoint[j] - positions[offset]);
				velocities[offset] = std::clamp(velocities[offset], -1.0, +1.0);
				double posNew = positions[offset] + velocities[offset];
				if (posNew < 0.0)
				{
					positions[offset] = 0.0;
					velocities[offset] = 0.0;
				}
				else if (posNew > 1.0)
				{
					positions[offset] = 1.0;
					velocities[offset] = 0.0;
				}
				else
				{
					positions[offset] = posNew;
				}
			}
			double value = modelCalc(positions + i * dims, dims);
			if (value < bestValues[i])
			{
				bestValues[i] = value;
				for (int j = 0; j < dims; j++)
				{
					int offset = i * dims + j;
					bestPositions[offset] = positions[offset];
				}
				if (value < globalBest)
				{
					globalBest = bestValues[i];
					std::cout << "new improvement to: " << globalBest << std::endl;
					for (int j = 0; j < dims; j++)
					{
						int offset = i * dims + j;
						globalBestPoint[j] = positions[offset];
						std::cout << globalBestPoint[j] << ", ";
					}
					std::cout << std::endl;
				}
			}
		}
		for (int j = 0; j < dims; j++)
		{
			int offset = (rand() % particles) * dims + j;
			probe[j] = bestPositions[offset];
		}
		double hybridValue = modelCalc(probe, dims);
		if (hybridValue < globalBest)
		{
			globalBest = hybridValue;
			std::cout << "new improvement to: " << globalBest << std::endl;
			for (int j = 0; j < dims; j++)
			{
				globalBestPoint[j] = probe[j];
				std::cout << globalBestPoint[j] << ", ";
			}
			std::cout << std::endl;
		}
		else
		{
			double worstValue = bestValues[0];
			int worst = 0;
			for (int j = 1; j < particles; j++)
			{
				if (worstValue < bestValues[j])
				{
					worstValue = bestValues[j];
					worst = j;
				}
			}
			if (hybridValue < worstValue)
			{
				bestValues[worst] = hybridValue;
				for (int j = 0; j < dims; ++j) {
					int offset = worst * dims + j;
					bestPositions[offset] = probe[j];
				}
			}
		}
	}
	delete[] positions;
	delete[] bestPositions;
	delete[] velocities;
	delete[] bestValues;
	delete[] probe;
	std::cout << "final error: " << globalBest << std::endl;
	for (int i = 0; i < dims; ++i) {
		std::cout << globalBestPoint[i] << ", ";
	}
	std::cout << std::endl;
	delete[] globalBestPoint;
}

int main(int argc, char* argv[]) {
	math::Matrix covar;
	covar.Load("D:\\Treasure\\Cpp\\ImageExperiments\\Data\\covar65x65.mat");
	const int N = static_cast<int>(sqrt(covar.Rows()));
	if (covar.Rows() != N * N || covar.Columns() != N * N) {
		std::cout << "expected square covariance matrix" << std::endl;
		return -1;
	}

	auto modelCalc = [&](double* modelParams, int dims) {
		return modelError(covar, N, modelParams);
		};
	particleSwarm(6, 20, 100, modelCalc);

	return 0;
}
