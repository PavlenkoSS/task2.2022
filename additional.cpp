#include "matrlib.h"

using namespace std;

double pol_1(double x)
{
	return x;
}
double abs_x(double x)
{
	return abs(x);
}

double runge_x(double x)
{
	return 1 / (25 * x * x + 1);
}
void fill_mat(string filename, matrix& x, matrix& f)
{
	ifstream fin;
	fin.open(filename);
	int i = 0,n=0;
	n = x.Dim_rows();
	fin >> n;
	double buf;
	for (i = 0; i < n; i++)
	{
		fin >> buf;
		x.SetMij(i, 0, buf);
	}
	for (i = 0; i < n; i++)
	{
		fin >> buf;
		f.SetMij(i, 0, buf);
	}
}
void uniform_points(double a, double b, double N, matrix& p)
{
	N = N - 1;
	double I = 0;
	for (int i = 0; i < N + 1; i++)
	{
		p.SetMij(i, 0, a + (b - a) * I / N);
		I++;
	}
}

void cheb_points(double a, double b, double N, matrix& p)
{
	double PI = 3.141592653589793;
	N = N + 1;
	double I = N;
	for (int i = N-1; i >= 0; i = i - 1)
	{
		p.SetMij(i, 0, (a + b) / 2 + (b - a) * cos(PI * (2 * I - 1) / (2 * N)) / 2);
		I--;
	}
}

// ����������� �����
void c_rand_points(double a, double b, double N, matrix& p)
{
	double I = 0;
	for (int i = 0; i < N; i++)
	{
		p.SetMij(i, 0, (double)(rand()) / RAND_MAX * (b - a) + a);
		I++;
	}
}	

void generation(string filename, double a, double b, double N, int p_type, double(*f)(double))
{
	ofstream fout;
	fout.open(filename);
	fout << N << endl;
	matrix p(N,1);
	switch (p_type)
	{
		case(2):
			cheb_points(a, b, N, p);
			break;
		case(3):
			c_rand_points(a, b, N, p);
			break;
		default:
			uniform_points(a, b, N, p);
			break;
	}
	for (int i = 0; i < N; i++)
	{
		fout << p.GetMij(i,0) << ' ';
	}
	fout << endl;
	for (int i = 0; i < N; i++)
	{
		fout << f(p.GetMij(i, 0)) << ' ';
	}
	fout.close();

}

double lagr_pol(double x, matrix& p, matrix& f) // �������
{
	double S = 0;
	double s = 1;
	int n = p.Dim_rows();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				s = s * (x - p.GetMij(j, 1)) / (p.GetMij(i, 1) - p.GetMij(j, 1));
			}	
		}
		S = S+f.GetMij(i, 1) * s;
		s = 1;
	}
	return S;
}

double stand_pol(double x, matrix& coef)
{
	double s = 0;
	int n = coef.Dim_rows();
	double X = 1;
	for (int i = 0; i < n+1; i++)
	{
		s = s + X * coef.GetMij(i, 0);
		X = X * x;
		//cout << i << ' ' << coef.GetMij(i, 1) << ' ' << X << ' ' << s <<'\n';
	}

	return s;
}

void val_filler(string filename, double a, double b, double(*f)(double))
{
	ofstream fout;
	fout.open(filename);
	for (double i = a; i < b; i += 1e-3)
	{
		fout << i << ' ' << f(i) << endl;
	} 
	fout.close();
}

void val_filler(string filename, double a, double b, double(*f)(double, matrix&), matrix& p)
{
	ofstream fout;
	fout.open(filename);
	for (double i = a; i < b; i += 1e-3)
	{
		fout << i << ' ' << f(i, p) << endl;
	}
	fout.close();
}

void val_filler(string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val)
{
	ofstream fout;
	fout.open(filename);
	for (double i = a; i < b; i += 1e-3)
	{
		fout << i << ' ' << f(i, p, val) << endl;
	}
	fout.close();
}	


void val_filler(int N, string filename, double a, double b, double(*f)(double))
{
	N = N - 1;
	ofstream fout;
	fout.open(filename);
	int I=0;
	for (double i = a; i < b; i +=  a + (b - a) * I / N)
	{
		fout << i << ' ' << f(i) << endl;
		I++;
	} 
	fout.close();
}

void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&), matrix& p)
{
	N = N - 1;
	ofstream fout;
	fout.open(filename);
	int I = 0;
	for (double i = a; i < b; i += a + (b - a) * I / N)
	{
		I++;
		fout << i << ' ' << f(i, p) << endl;
	}
	fout.close();
}

void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val)
{
	N = N - 1;
	ofstream fout;
	fout.open(filename);
	int I = 0;
	for (double i = a; i < b; i +=  a + (b - a) * I / N)
	{
		I++;
		fout << i << ' ' << f(i, p, val) << endl;
	}
	fout.close();
}	
