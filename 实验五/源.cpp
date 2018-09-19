
#include <iostream>
#include <math.h>

using namespace::std;

double** gen_h(int n) {
	double** H = new double*[n];
	for(int i = 0; i < n; i++)
		H[i] = new double[n];
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			H[i][j] = 1.0 / (i + j + 1);
	return H;
}

double* mul(double** A, double* x,int n) {
	double* b = new double[n];
	for(int i = 0; i < n; i++) {
		b[i] = 0.0;
		for(int j = 0; j < n; j++) {
			b[i] += A[i][j] * x[j];
		}
	}
	return b;
}

double** trans(double** A, int n) {
	double** H = new double*[n];
	for(int i = 0; i < n; i++)
		H[i] = new double[n];
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			H[i][j] = A[j][i];
	return H;
}

double norm_1(double* x, int n) {
	double result = 0.0;
	for(int i = 0; i < n; i++) {
		if(result < x[i])
			result = x[i];
	}
	return result;
}


double norm_1(double** A,int n) {
	double result = 0.0;
	for(int i = 0; i < n; i++) {
		double sum = 0.0;
		for(int j = 0; j < n; j++)
			sum += abs(A[i][j]);
		if(result < sum)
			result = sum;
	}
	return result;
}

void LU_Descomposition(double** A, double** L, double** U, int n) {
	for(int i = 0; i < n; i++) {
		U[0][i] = A[0][i];
		L[i][0] = A[i][0] / U[0][0];
	}
	for(int r = 1; r < n; r++) {
		for(int i = r; i < n; i++) {
			double s = 0.0, t = 0.0;
			for(int k = 0; k < r; k++) {
				s += L[r][k] * U[k][i];
				t += L[i][k] * U[k][r];
			}
			U[r][i] = A[r][i] - s;
			L[i][r] = (A[i][r] - t) / U[r][r];
		}
		L[r][r] = 1.0;
	}
}

double** inverse(double** A,int n) {
	double** L = new double*[n];
	for(int i = 0; i < n; i++)
		L[i] = new double[n];
	double** U = new double*[n];
	for(int i = 0; i < n; i++)
		U[i] = new double[n];
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			L[i][j] = U[i][j] = 0.0;
	LU_Descomposition(A, L, U, n);
	double** B = new double*[n];
	double** Y = new double*[n];
	double** X = new double*[n];
	for(int j = 0; j < n; j++) {
		B[j] = new double[n];
		Y[j] = new double[n];
		X[j] = new double[n];
		double *b = B[j], *y = Y[j], *x = X[j];
		for(int i = 0; i < n; i++)
			b[i] = 0.0;
		b[j] = 1.0;
		y[0] = b[0];
		for(int i = 1; i < n; i++) {
			double t = 0.0;
			for(int k = 0; k < i; k++)
				t += L[i][k] * y[k];
			y[i] = b[i] - t;
		}
		x[n - 1] = y[n - 1] / U[n - 1][n - 1];
		for(int i = n - 2; i >= 0; i--) {
			double s = 0.0;
			for(int k = i + 1; k < n; k++)
				s += U[i][k] * x[k];
			x[i] = (y[i] - s) / U[i][i];
		}
	}
	double** _X = trans(X, n);
	return _X;
}

double cond(double** A, int n,double(*norm)(double**, int)) {
	return norm(A, n)*norm(inverse(A,n), n);
}

double task2(int n) {
	double** H = gen_h(n);
	double* x = new double[n];
	for(int i = 0; i < n; i++)
		x[i] = 1.0;
	double* b = mul(H, x, n);
	double** L = new double*[n];
	for(int i = 0; i < n; i++)
		L[i] = new double[n];
	double** U = new double*[n];
	for(int i = 0; i < n; i++)
		U[i] = new double[n];
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			L[i][j] = U[i][j] = 0.0;
	LU_Descomposition(H, L, U, n);
	double* y = new double[n];
	double* _x = new double[n];
	y[0] = b[0];
	for(int i = 1; i < n; i++) {
		double t = 0.0;
		for(int k = 0; k < i; k++)
			t += L[i][k] * y[k];
		y[i] = b[i] - t;
	}
	_x[n - 1] = y[n - 1] / U[n - 1][n - 1];
	for(int i = n - 2; i >= 0; i--) {
		double s = 0.0;
		for(int k = i + 1; k < n; k++)
			s += U[i][k] * _x[k];
		_x[i] = (y[i] - s) / U[i][i];
	}
	//for(int i = 0; i < n; i++)
	//	cout << _x[i] << endl;
	double* r = new double[n];
	double* Hn_x = mul(H, _x, n);
	for(int i = 0; i < n; i++)
		r[i] = b[i] - Hn_x[i];
	cout << "norm_1(r): " << norm_1(r, n) << endl;
	double* dx = new double[n];
	for(int i = 0; i < n; i++)
		dx[i] = _x[i] - x[i];
	cout << "norm_1(dx): " << norm_1(dx, n) << endl;
	return  norm_1(dx, n);
}

void task4() {
	int n = 10;
	double** H = gen_h(n);
	double* x = new double[n];
	for(int i = 0; i < n; i++)
		x[i] = 1.0;
	double* b = mul(H, x, n);
	b[0] += 0.0000001;
	double** L = new double*[n];
	for(int i = 0; i < n; i++)
		L[i] = new double[n];
	double** U = new double*[n];
	for(int i = 0; i < n; i++)
		U[i] = new double[n];
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			L[i][j] = U[i][j] = 0.0;
	LU_Descomposition(H, L, U, n);
	double* y = new double[n];
	double* _x = new double[n];
	y[0] = b[0];
	for(int i = 1; i < n; i++) {
		double t = 0.0;
		for(int k = 0; k < i; k++)
			t += L[i][k] * y[k];
		y[i] = b[i] - t;
	}
	_x[n - 1] = y[n - 1] / U[n - 1][n - 1];
	for(int i = n - 2; i >= 0; i--) {
		double s = 0.0;
		for(int k = i + 1; k < n; k++)
			s += U[i][k] * _x[k];
		_x[i] = (y[i] - s) / U[i][i];
	}
	double* r = new double[n];
	double* Hn_x = mul(H, _x, n);
	for(int i = 0; i < n; i++)
		r[i] = b[i] - Hn_x[i];
	cout << "norm_1(r): " << norm_1(r, n) << endl;
	double* dx = new double[n];
	for(int i = 0; i < n; i++)
		dx[i] = _x[i] - x[i];
	cout << "norm_1(dx): " << norm_1(dx, n) << endl;
}

void task5() {
	int n = 5;
	double dx = 0.0;
	while(dx < 1.0) {
		cout << n << endl;
		dx = task2(n++);
	}
}

int main() {
	cout << "cond(H3): " << cond(gen_h(3), 3, norm_1) << endl;
	cout << "cond(H4): " << cond(gen_h(4), 4, norm_1) << endl;
	task2(10);
	task4();
	task5();
	int n;
	cin >> n;
}