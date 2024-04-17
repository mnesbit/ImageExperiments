#include "../inc/mathvector.h"
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <fstream>
#include "../inc/mathmatrix.h"

namespace math
{

template <class T>
inline T abs(const T &x) { return (x > 0 ? x : -x); };

template <class T>
inline int sign(const T &x) { return (x >= 0 ? 1 : -1); };

template <class T>
inline T square(const T &x) { return x*x; };

Vector::Vector()
	: m(0)
	, data(nullptr)
{
}

Vector::Vector(size_t length,const double* vectordata)
	: m(length)
	, data(new double[length])
{
	if(vectordata)
	{
		memcpy(data,vectordata,length*sizeof(double));
	}
	else
	{
		memset(data,0,length*sizeof(double));
	}
}

Vector::Vector(const Vector& vector)
	: m(vector.Length())
	, data(new double[m])
{
	memcpy(data, vector.data, m * sizeof(double));
}

Vector::Vector(Vector&& vector) noexcept
	: m(vector.m)
	, data(vector.data)
{
	vector.m = 0;
	vector.data = nullptr;
}

Vector::~Vector()
{
	delete[] data;
	data = nullptr;
}

Vector& Vector::operator=(const Vector& vector)
{
	if (this != &vector)
	{
		delete[] data;
		m = vector.m;
		data = new double[m];
		memcpy(data, vector.data, m * sizeof(double));
	}
	return *this;
}

Vector& Vector::operator=(Vector&& vector) noexcept
{
	if (this != &vector)
	{
		delete[] data;
		m = vector.m;
		data = vector.data;
		vector.m = 0;
		vector.data = nullptr;
	}
	return *this;
}

double& Vector::operator[](size_t nSubscript)
{
	if((nSubscript>=0) && (nSubscript<m))
	{
		return data[nSubscript];
	}
	throw new std::range_error("invalid vector index");
}

const double& Vector::operator[](size_t nSubscript) const
{
	if((nSubscript>=0) && (nSubscript<m))
	{
		return data[nSubscript];
	}
	throw new std::range_error("invalid vector index");
}

void Vector::Add(const Vector& right)
{
	if(Length()!=right.Length())throw  new std::range_error("mismatched sizes");
	double* pDataL=data;
	double* pDataR=right.data;
	for(size_t i=0;i<m;i++)
	{
		double value=(*pDataL)+(*pDataR++);
		*pDataL++=value;
	}
}

void Vector::Subtract(const Vector& right)
{
	if(Length()!=right.Length())throw  new std::range_error("mismatched sizes");
	double* pDataL=data;
	double* pDataR=right.data;
	for(size_t i=0;i<m;i++)
	{
		double value=(*pDataL)-(*pDataR++);
		*pDataL++=value;
	}
}

void Vector::ReverseSubtract(const Vector& left)
{
	if(Length()!=left.Length())throw  new std::range_error("mismatched sizes");
	double* pDataR=data;
	double* pDataL=left.data;
	for(size_t i=0;i<m;i++)
	{
		double value=(*pDataL++)-(*pDataR);
		*pDataR++=value;
	}
}

void Vector::Scale(double scaling)
{
	double* pData=data;
	for(size_t i=0;i<m;i++)
	{
		double value=scaling*(*pData);
		*pData++=value;
	}
}

void Vector::AddConstant(double constant)
{
	double* pData=data;
	for(size_t i=0;i<m;i++)
	{
		double value=constant+(*pData);
		*pData++=value;
	}
}

std::ostream& operator<<(std::ostream& os, const Vector& vector)
{
	const double *pData=vector.Data();
	for(size_t i=0;i<vector.Length();i++)
	{
		os<<(*pData++)<<std::endl;
	}
	os <<std::endl;
	return os;
}


double Vector::Magnitude()
{
	double* pData=data;
	double biggest=0.0; //reduce overflow by diving through
	for(size_t i=0;i<m;i++)
	{
		if(abs(*pData)>biggest)
		{
			biggest=abs(*pData);
		}
		pData++;
	}
	pData=data;
	double tot=0.0;
	for(size_t i=0;i<m;i++)
	{
		tot+=square((*pData)/biggest);
		pData++;
	}
	return biggest*sqrt(tot);
}

double Vector::SquaredMagnitude()
{
	double tot=0.0;
	double* pData=data;
	for(size_t i=0;i<m;i++)
	{
		tot+=square(*pData);
		pData++;
	}
	return tot;
}

double Vector::DotProduct(const Vector& right)
{
	if(Length()!=right.Length())throw  new std::range_error("mismatched sizes");
	double tot=0.0;
	double* pDataL=data;
	double* pDataR=right.data;
	for(size_t i=0;i<m;i++)
	{
		tot+=(*pDataL++)*(*pDataR++);
	}
	return tot;
}


void Vector::Save(const char* filename) const
{
	std::ofstream file(filename, std::ios::out | std::ios::binary);
	if(file)
	{
		file << "VEC ";
		file.write((char*)&m,sizeof(m));
		file.write((char*)data,m*sizeof(double));
	}
	else
	{
		throw new std::runtime_error("IO error");
	}
}

void Vector::Load(const char* filename)
{
	std::ifstream file(filename,std::ios::in|std::ios::binary);
	if(file)
	{
		char sigbuff[5];
		file.read(sigbuff,4);
		sigbuff[4]=0;
		if(strcmp(sigbuff,"VEC "))
		{
			throw new std::runtime_error("IO error");
		}
		delete data;
		data=0;
		file.read((char*)&m,sizeof(m));
		data=new double[m];
		file.read((char*)data,m*sizeof(double));
	}
	else
	{
		throw new std::runtime_error("IO error");
	}
}

void Vector::Reset()
{
	delete [] data;
	data=0;
	m=0;
}

Vector Add(const Vector& left,const Vector& right)
{
	if((left.Length()!=right.Length()))throw new std::range_error("mismatched sizes");
	Vector result(left);
	result.Add(right);
	return result;
}

Vector Subtract(const Vector& left,const Vector& right)
{
	if((left.Length()!=right.Length()))throw new std::range_error("mismatched sizes");
	Vector result(left);
	result.Subtract(right);
	return result;
}


}
