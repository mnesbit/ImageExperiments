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
	double& operator[](int nSubscript);
	const double& operator[](int nSubscript) const;
	int Length() const { return m; }
	~MatrixRow();
private:
	MatrixRow(int rowlength,double* rowdata); //only constructable by Matrix
	MatrixRow& operator=(const MatrixRow& matrix); // not assignable
	int m;
	double* rowdata;
};

class Matrix
{
public:
	friend class Vector;
	friend class MatrixRow;
	Matrix();
	Matrix(int rows,int columns,const double* arraydata=0);
	Matrix(int squarewidth,const double* arraydata=0);
	Matrix(const Matrix& matrix);
	Matrix(const Vector& vector);
	~Matrix();

	Matrix& operator=(const Matrix& matrix);
	MatrixRow operator[](int nSubscript);
	const MatrixRow operator[](int nSubscript) const;

	int Rows() const { return m; }
	int Columns() const { return n; }
	Vector GetDiag() const;
	void SetDiag(const Vector& diag);
	Vector GetRow(int nSubscript) const;
	Vector GetColumn(int nSubscript) const;
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
	int m;
	int n;
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
