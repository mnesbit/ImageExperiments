#include "../inc/mathmatrix.h"
#include <cstring>
#include <stdexcept>
#include <fstream>
#include "../inc/mathvector.h"
#include <algorithm>

namespace math
{

	
template <class T>
inline T abs(const T &x) { return (x > 0 ? x : -x); };

template <class T>
inline int sign(const T &x) { return (x >= 0 ? 1 : -1); };

template <class T>
inline T square(const T &x) { return x*x; };

Matrix::Matrix()
	: m(0)
	, n(0)
	, data(0)
{
}

Matrix::Matrix(size_t rows, size_t columns,const double* arraydata)
	: m(rows)
	, n(columns)
	, data(new double[m*n])
{
	if(arraydata)
	{
		memcpy(data,arraydata,m*n*sizeof(double));
	}
	else
	{
		memset(data,0,m*n*sizeof(double));
	}
}

Matrix::Matrix(size_t squarewidth,const double* arraydata)
	: m(squarewidth)
	, n(squarewidth)
	, data(new double[m*n])
{
	if(arraydata)
	{
		memcpy(data,arraydata,m*n*sizeof(double));
	}
	else
	{
		memset(data,0,m*n*sizeof(double));
	}
}

Matrix::Matrix(const Matrix& matrix)
	: m(matrix.Rows())
	, n(matrix.Columns())
	, data(new double[m * n])
{
	memcpy(data, matrix.data, m * n * sizeof(double));
}

Matrix::Matrix(Matrix&& matrix) noexcept
	: m(0)
	, n(0)
	, data(nullptr)
{
	m = matrix.m;
	n = matrix.n;
	data = matrix.data;
	matrix.m = 0;
	matrix.n = 0;
	matrix.data = nullptr;
}

Matrix::Matrix(const Vector& vector)
	: m(vector.Length())
	, n(vector.Length())
	, data(new double[m*n])
{
	for(size_t i=0;i<m;i++)
	{
		data[i*(n+1)]=vector[i];
	}
}

Matrix::~Matrix()
{
	delete[] data;
	data = nullptr;
}

Matrix& Matrix::operator=(const Matrix& matrix)
{
	if (this != &matrix)
	{
		delete[] data;
		m = matrix.m;
		n = matrix.n;
		data = new double[m * n];
		memcpy(data, matrix.data, m * n * sizeof(double));
	}
	return *this;
}

Matrix& Matrix::operator=(Matrix&& matrix) noexcept
{
	if (this != &matrix)
	{
		delete[] data;
		m = matrix.m;
		n = matrix.n;
		data = matrix.data;
		matrix.m = 0;
		matrix.n = 0;
		matrix.data = nullptr;
	}
	return *this;
}

Vector Matrix::GetDiag() const
{
	Vector diag(std::min(m,n));
	for(size_t i=0;i<std::min(m,n);i++)
	{
		diag[i]=data[i*(n+1)];
	}
	return diag;
}

void Matrix::SetDiag(const Vector& diag)
{
	if(std::min(m,n)!=diag.Length())throw new std::range_error("invalid vector length");
	for(size_t i=0;i<std::min(m,n);i++)
	{
		data[i*(n+1)]=diag[i];
	}
}

Vector Matrix::GetRow(size_t nSubscript) const
{
	if((nSubscript>=0)&&(nSubscript<m))
	{
		Vector row(n);
		double *pDataI=&data[n*nSubscript];
		double *pDataO=row.Data();
		for(size_t i=0;i<n;i++)
		{
			*pDataO=*pDataI;
			pDataO++;
			pDataI++;
		}
		return row;
	}
	throw new std::range_error("invalid row index");
}

Vector Matrix::GetColumn(size_t nSubscript) const
{
	if((nSubscript>=0)&&(nSubscript<n))
	{
		Vector column(m);
		double *pDataI=&data[nSubscript];
		double *pDataO=column.Data();
		for(size_t i=0;i<m;i++)
		{
			*pDataO=*pDataI;
			pDataO++;
			pDataI+=n;
		}
		return column;
	}
	throw new std::range_error("invalid column index");
}

Matrix Matrix::GetTranspose() const
{
	Matrix transp(n,m);
	double* tdata=transp.Data();
	for(size_t i=0;i<m;i++)
	{
		for(size_t j=0;j<n;j++)
		{
			tdata[j*m+i]=data[i*n+j];
		}
	}
	return transp;
}

MatrixRow Matrix::operator[](size_t nSubscript)
{
	if((nSubscript>=0)&&(nSubscript<m))
	{
		return MatrixRow(n,&data[nSubscript*n]);
	}
	throw new std::range_error("invalid row index");
}

const MatrixRow Matrix::operator[](size_t nSubscript) const
{
	if((nSubscript>=0)&&(nSubscript<m))
	{
		return MatrixRow(n,&data[nSubscript*n]);
	}
	throw new std::range_error("invalid row index");
}

MatrixRow::MatrixRow(size_t rowlength,double* rowdata)
	: m(rowlength)
	, rowdata(rowdata)
{
}

MatrixRow::~MatrixRow()
{
}

double& MatrixRow::operator[](size_t nSubscript)
{
	if((nSubscript>=0) && (nSubscript<m))
	{
		return rowdata[nSubscript];
	}
	throw new std::range_error("invalid row index");
}

const double& MatrixRow::operator[](size_t nSubscript) const
{
	if((nSubscript>=0) && (nSubscript<m))
	{
		return rowdata[nSubscript];
	}
	throw new std::range_error("invalid row index");
}

void Matrix::MakeIdentityMatrix()
{
	if(m!=n)throw new std::range_error("matrix must be square");
	memset(data,0,m*n*sizeof(double));
	for(size_t i=0;i<m;i++)
	{
		data[i*(n+1)]=1.0;
	}
}

void Matrix::Add(const Matrix& right)
{
	if((Rows()!=right.Rows())||(Columns()!=right.Columns()))throw new std::range_error("mismatched sizes");
	double* pDataL=data;
	double* pDataR=right.data;
	for(size_t i=0;i<(m*n);i++)
	{
		double value=(*pDataL)+(*pDataR++);
		*pDataL++=value;
	}
}

void Matrix::Subtract(const Matrix& right)
{
	if((Rows()!=right.Rows())||(Columns()!=right.Columns()))throw new std::range_error("mismatched sizes");
	double* pDataL=data;
	double* pDataR=right.data;
	for(size_t i=0;i<(m*n);i++)
	{
		double value=(*pDataL)-(*pDataR++);
		*pDataL++=value;
	}
}

void Matrix::ReverseSubtract(const Matrix& left)
{
	if((Rows()!=left.Rows())||(Columns()!=left.Columns()))throw new std::range_error("mismatched sizes");
	double* pDataR=data;
	double* pDataL=left.data;
	for(size_t i=0;i<(m*n);i++)
	{
		double value=(*pDataR)-(*pDataL++);
		*pDataR++=value;
	}
}

void Matrix::Scale(double scaling)
{
	double* pData=data;
	for(size_t i=0;i<(m*n);i++)
	{
		double value=scaling*(*pData);
		*pData++=value;
	}
}

void Matrix::AddConstant(double constant)
{
	double* pData=data;
	for(size_t i=0;i<(m*n);i++)
	{
		double value=constant+(*pData);
		*pData++=value;
	}
}

Matrix Matrix::GetSquare()const
{
	Matrix result(m,m);
	double* pDataO=result.Data();
	for(size_t i=0;i<m;i++)
	{
		double tot=0.0;
		double* pDataL=&data[i*n];
		double* pDataR;
		for(size_t k=0;k<n;k++)
		{
			tot+=square(*pDataL);
			pDataL++;
		}
		pDataO[i*(m+1)]=tot;
		for(size_t j=0;j<i;j++)
		{
			tot=0.0;
			pDataL=&data[i*n];
			pDataR=&data[j*n];
			for(size_t k=0;k<n;k++)
			{
				tot+=(*pDataL)*(*pDataR);
				pDataL++;
				pDataR++;
			}
			pDataO[i*m+j]=tot;
			pDataO[j*m+i]=tot;
		}
	}
	return result;
}

void Matrix::Save(const char* filename) const
{
	std::ofstream file(filename, std::ios::out | std::ios::binary);
	if(file)
	{
		file << "MTX ";
		file.write((char*)&m,sizeof(m));
		file.write((char*)&n,sizeof(n));
		file.write((char*)data,m*n*sizeof(double));
	}
	else
	{
		throw new std::runtime_error("IO error");
	}
}

void Matrix::Load(const char* filename)
{
	std::ifstream file(filename,std::ios::in|std::ios::binary);
	if(file)
	{
		char sigbuff[5];
		file.read(sigbuff,4);
		sigbuff[4]=0;
		if(strcmp(sigbuff,"MTX "))
		{
			throw new std::runtime_error("IO error");
		}
		delete data;
		data=0;
		file.read((char*)&m,sizeof(m));
		file.read((char*)&n,sizeof(n));
		data=new double[m*n];
		file.read((char*)data,m*n*sizeof(double));
	}
	else
	{
		throw new std::runtime_error("IO error");
	}
}


Matrix Add(const Matrix& left,const Matrix& right)
{
	if((left.Rows()!=right.Rows())||(left.Columns()!=right.Columns()))throw new std::range_error("mismatched sizes");
	Matrix result(left);
	result.Add(right);
	return result;
}

Matrix Subtract(const Matrix& left,const Matrix& right)
{
	if((left.Rows()!=right.Rows())||(left.Columns()!=right.Columns()))throw new std::range_error("mismatched sizes");
	Matrix result(left);
	result.Subtract(right);
	return result;
}

Matrix Multiply(const Matrix& left,const Matrix& right)
{
	if(left.Columns()!=right.Rows())throw  new std::range_error("mismatched sizes");
	size_t nrows=left.Rows();
	size_t ncols=left.Columns();
	size_t mcols=right.Columns();
	Matrix result(nrows,mcols);
	const double* pDataL=left.Data();
	const double* pDataR=right.Data();
	double* pDataO=result.Data();
	const double *pL,*pR;
	for(size_t i=0;i<nrows;i++)
	{
		for(size_t j=0;j<mcols;j++)
		{
			double tot=0.0;
			pL=pDataL;
			pR=pDataR;
			for(size_t k=0;k<ncols;k++)
			{
				tot+=(*pL)*(*pR);
				pL++;
				pR+=mcols;
			}
			*pDataO++=tot;
			pDataR++;
		}
		pDataL+=ncols;
		pDataR=right.Data();
	}
	return result;
}

Vector Multiply(const Matrix& left,const Vector& right)
{
	if(left.Columns()!=right.Length())throw  new std::range_error("mismatched sizes");
	Vector result(left.Rows());
	const double *pDataL=left.Data();
	const double *pDataR=right.Data();
	double *pDataO=result.Data();
	const double *pL,*pR;
	for(size_t i=0;i<left.Rows();i++)
	{
		double tot=0.0;
		pL=pDataL;
		pR=pDataR;
		for(size_t j=0;j<right.Length();j++)
		{
			tot+=(*pL)*(*pR);
			pL++;
			pR++;
		}
		*pDataO++=tot;
		pDataL+=left.Columns();
	}
	return result;
}

void Matrix::Reset()
{
	delete[] data;
	data=0;
	m=0;
	n=0;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix)
{
	const double *pData=matrix.Data();
	for(size_t i=0;i<matrix.Rows();i++)
	{
		for(size_t j=0;j<matrix.Columns();j++)
		{
			os<<(*pData++)<<" ";
		}
		os <<std::endl;
	}
	os <<std::endl;
	return os;
}


}

