#ifndef _MATH_VECTOR_H_
#define _MATH_VECTOR_H_

#include "framework.h"

#include<ostream>

namespace math
{

class Matrix;

class Vector
{
public:
	friend class Matrix;
	Vector();
	Vector(size_t length,const double* vectordata=0);
	Vector(const Vector& vector);
	Vector(Vector&& vector) noexcept;
	~Vector();

	Vector& operator=(const Vector& vector);
	Vector& operator=(Vector&& vector) noexcept;
	const double& operator[](size_t nSubscript) const;
	double& operator[](size_t nSubscript);

	size_t Length() const { return m; }
	double* Data() { return data; }
	const double* Data() const { return data; }

	void Add(const Vector& right);
	void Subtract(const Vector& right);
	void ReverseSubtract(const Vector& left);
	void Scale(double scaling);
	void AddConstant(double constant);
	double Magnitude();
	double SquaredMagnitude();
	double DotProduct(const Vector& right);
	void Save(const char* filename) const;
	void Load(const char* filename);
	void Reset();
private:
	size_t m;
	double* data;
};

std::ostream& operator<<(std::ostream& os, const Vector& vector);
Vector Add(const Vector& left,const Vector& right);
Vector Subtract(const Vector& left,const Vector& right);

}

#endif //_MATH_VECTOR_H_
