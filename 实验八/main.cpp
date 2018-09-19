#include <iostream>
#include <math.h>
#include <iomanip>
using namespace::std;

double eps = 0.00001;

double max(int dim, double* l) {
	double tem = l[0];
	for(int i = 0; i < dim; i++)
		if(l[i] > tem)
			tem = l[i];
	return tem;
}

double* multiply(int dim, double* u, double** A) {
	double* result = new double[dim];
	for(int i = 0; i < dim; i++) {
		double tem = 0.0;
		for(int j = 0; j < dim; j++)
			tem += A[i][j] * u[j];
		result[i] = tem;
	}
	return result;
}

double* divide(int dim, double* v, double m) {
	double* result = new double[dim];
	for(int i = 0; i < dim; i++)
		result[i] = v[i] / m;
	return result;
}

void power_method(int dim, double** A) {
	double* v = new double[dim];
	double* u = new double[dim];
	for(int i = 0; i < dim; i++)
		v[i] = u[i] = 1.0;
	double _m, m = 0.0;
	do {
		_m = m;
		v = multiply(dim, u, A);
		m = max(dim, v);
		u = divide(dim, v, m);
	} while(fabs(m - _m) > eps);
	cout << "lamda1: " << m << endl;
	cout << "x1: " << endl;
	for(int i = 0; i < dim; i++)
		cout << u[i] << " ";
	cout << endl;
}
double _A[3][3] = { { 5,-4,1 }, { -4,6,-4 }, { 1,-4,7 } };
double _B[4][4] = { { 25,-41,10,-4 }, { -41,68,-17,10 }, { 10,-17,5,-3 }, { -6,10,-3,2 } };

int main() {
	double** A = new double*[3];
	for(int i = 0; i<3; i++)
		A[i] = new double[3];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			A[i][j] = _A[i][j];
	double** B = new double*[4];
	for(int i = 0; i<4; i++)
		B[i] = new double[4];
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			B[i][j] = _B[i][j];
	power_method(3, A);
	power_method(4, B);
	int n;
	cin >> n;
}