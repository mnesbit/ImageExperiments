#ifndef _SYMM_EIGEN_H_
#define _SYMM_EIGEN_H_

#include "framework.h"

#include "mathvector.h"
#include "mathmatrix.h"

namespace math
{

class Vector;
class Matrix;

class SymmetricEigenDecomposition
{
public:
	SymmetricEigenDecomposition(const Matrix& input);
	~SymmetricEigenDecomposition();
	Vector EigenValues()const { return m_eigenValues; }
	Matrix EigenVectors()const { return m_eigenVectors; }
private:
	Matrix m_eigenVectors;
	Vector m_eigenValues;
};

}

#endif //_SYMM_EIGEN_H_
