#ifndef matrlib
#define matrlib
#include <iostream>
#include <vector>
#include <fstream>
#include <typeinfo>
#include <ctime>
#include<math.h>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>

using namespace std;


class matrix
{
private:
	int n; // строк 
	int m; // столбцов
	double* mat;
public:
	matrix()
	{
		n = 0;
		m = 0;
		mat = nullptr;
	}
	matrix(int m_rows, int m_cols)
	{
		n = m_rows;
		m = m_cols;
		mat = new double[n * m];
	}
	matrix(int m_rows)
	{
		n = m_rows;
		m = m_rows;
		mat = new double[n * m];
	}
	int Dim_rows()
	{
		return n;
	}
	int Dim_cols()
	{
		return m;
	}
	double GetMij(int i, int j)
	{
		if ((m > 0) && (n > 0))
			return mat[i * m + j];
		else
			return 0;
	}
	void SetMij(int i, int j, double value)
	{
		if ((i < 0) || (i >= n))
			return;
		if ((j < 0) || (j >= m))
			return;
		mat[i * m + j] = value;
	}

	matrix operator=(const matrix& M)
	{
		if ((n > 0) && (m > 0))
		{
			delete[] mat;
		}

		m = M.m;
		n = M.n;
		mat = new double[n * m];

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				mat[i * n + j] = M.mat[i * n + j];
		return *this;
	}

	void out_mat()
	{

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				cout << setprecision(3);
				cout << mat[i * m + j] << ' ';
			}
			cout << endl;
		}
		cout << endl << endl;
	}
	int fil_mat_file(string filename)
	{
		ifstream fin;
		fin.open(filename);
		if (!fin.is_open())
		{
			cout << "File was not opened" << endl;
			return -1;
		}
		for (int i = 0; i < n * m ; i++)
		{
			if (!fin.eof())
			{
				fin >> mat[i];
			}
			/*if ((fin.eof()) && (i != (n * m - 1)))
			{
				cout << "Filling error" << endl;
				return -1;
			}
			if ((!fin.eof()) && (i == (n * m - 1)))
			{
				cout << "Filling error" << endl;
				return -1;
			}*/
		}
		fin.close();
		return 0;
	}
	void fil_mat_by_x(const matrix& x)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				mat[i * m + j] = pow(x.mat[i], j);
	}
	void fil_mat_by_exp(const matrix& x)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				mat[i * m + j] = exp(x.mat[i]);
	}
	void friend up_triangle(matrix& A, matrix& b)
	{
		double d;
		int n = A.Dim_rows();
		int m = A.Dim_cols();
		for (int k = 0; k < n; k++) // ������ ���
		{
			if (abs(A.mat[k * m + k]) > 1e-14)
			{
				for (int i = k + 1; i < m; i++)
				{
					if (abs(A.mat[i * m + k]) > 1e-14)
					{
						d = A.mat[i * m + k] / A.mat[k * m + k];
						for (int j = k; j < n; j++)
						{
							A.mat[i * m + j] = A.mat[i * m + j] - d * A.mat[k * m + j];
						}
						b.mat[i] = b.mat[i] - d * b.mat[k];
					}
				}
			}
			else
			{
				double buf = 0;
				for (int i = k + 1; i < m; i++)
				{
					if (abs(A.mat[i * m + k]) > 1e-14)
					{
						for (int j = k; j < m; j++)
						{
							buf = A.mat[k * m + j];
							A.mat[k * m + j] = A.mat[i * m + j];
							A.mat[i * m + j] = buf;
						}
					}
				}
				for (int i = k + 1; i < m+1; i++)
				{
					if (i == m + 1)
					{
						cout << "System can not be solved";
						break;
					}
					if (abs(A.mat[i * m + k]) > 1e-14)
					{
						d = A.mat[i * m + k] / A.mat[k * m + k];
						for (int j = k; j < n; j++)
						{
							A.mat[i * m + j] = A.mat[i * m + j] - d * A.mat[k * m + j];
						}
						b.mat[i] = b.mat[i] - d * b.mat[k];
					}
				}
			}
		}
	}
	void friend gauss_back(matrix& A, matrix& b)
	{
		double d;
		int n = A.Dim_rows();
		int m = A.Dim_cols();
		for (int i = 0; i < n; i++)
		{
			d = A.mat[i * m + i];
			for (int j = i; j < n; j++)
			{
				A.mat[i * m + j] /= d;
			}
			b.mat[i] /= d;
		}
		for (int j = n - 1; j >= 0; j--)
		{
			for (int i = 0; i < j; i++)
			{
				b.mat[i] = b.mat[i] - b.mat[j] * A.mat[i * m + j];
				A.mat[i * m + j] = 0;
			}
		}
		//x = b;
	}
	/*~matrix()
	{
		delete[] mat;
	}*/

};
void cheb_points(double a, double b, double N, matrix& p);
double stand_pol(double x, matrix& coef);
double lagr_pol(double x, matrix& p, matrix& f); // �������
void uniform_points(double a, double b, double N, matrix& p);
void c_rand_points(double a, double b, double N, matrix& p);
void generation(string filename, double a, double b, double N, int p_type, double(*f)(double));
double abs_x(double x);
double runge_x(double x); void fill_mat(string filename, matrix& x, matrix& f);
void val_filler(string filename, double a, double b, double(*f)(double));
void val_filler(string filename, double a, double b, double(*f)(double, matrix&), matrix& p);
void val_filler(string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val);
double pol_1(double x);
void val_filler(int N, string filename, double a, double b, double(*f)(double));
void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&), matrix& p);
void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val);
#endif
