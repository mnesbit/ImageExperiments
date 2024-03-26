#ifndef _MATH_MATRIX_H_
#define _MATH_MATRIX_H_

#include "framework.h"

#include <ostream>

namespace math
{

class Vector;

class MatrixRow
{
public:
	friend class Matrix;
	double& operator[](size_t nSubscript);
	const double& operator[](size_t nSubscript) const;
	size_t Length() const { return m; }
	~MatrixRow();
private:
	MatrixRow(size_t rowlength,double* rowdata); //only constructable by Matrix
	MatrixRow& operator=(const MatrixRow& matrix); // not assignable
	size_t m;
	double* rowdata;
};

class Matrix
{
public:
	friend class Vector;
	friend class MatrixRow;
	Matrix();
	Matrix(size_t rows, size_t columns,const double* arraydata=0);
	Matrix(size_t squarewidth,const double* arraydata=0);
	Matrix(const Matrix& matrix);
	Matrix(Matrix&& matrix) noexcept;
	Matrix(const Vector& vector);
	~Matrix();

	Matrix& operator=(const Matrix& matrix);
	Matrix& operator=(Matrix&& matrix) noexcept;
	MatrixRow operator[](size_t nSubscript);
	const MatrixRow operator[](size_t nSubscript) const;

	size_t Rows() const { return m; }
	size_t Columns() const { return n; }
	Vector GetDiag() const;
	void SetDiag(const Vector& diag);
	Vector GetRow(size_t nSubscript) const;
	Vector GetColumn(size_t nSubscript) const;
	Matrix GetTranspose() const;
	double* Data() { return data; }
	const double* Data() const { return data; }
	void MakeIdentityMatrix();
	void Add(const Matrix& right);
	void Subtract(const Matrix& right);
	void ReverseSubtract(const Matrix& left);
	void Scale(double scaling);
	void AddConstant(double constant);
	Matrix GetSquare()const;
	void Save(const char* filename) const;
	void Load(const char* filename);
	void Reset();
private:
	size_t m;
	size_t n;
	double* data;
};
std::ostream& operator<<(std::ostream& os, const Matrix& matrix);

Matrix Add(const Matrix& left,const Matrix& right);
Matrix Subtract(const Matrix& left,const Matrix& right);
Matrix Multiply(const Matrix& left,const Matrix& right);
Matrix Multiply(const Vector& left,const Vector& right);
Vector Multiply(const Matrix& left,const Vector& right);

}

#endif //_MATH_MATRIX_H_
