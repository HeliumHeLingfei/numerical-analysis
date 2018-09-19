#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

int m = 7;
double x[7] = { -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0 };
double y[7] = { -4.467,-0.452, 0.551, 0.048, -0.447, 0.549, 4.552 };

double Alpha(double* i) {
	double temp1 = 0.0, temp2 = 0.0, temp = 0.0;;
	for(int k = 0; k < m; k++) {
		temp = i[k] * i[k];
		temp1 += temp * x[k];
		temp2 += temp;
	}
	return temp1/temp2;
}

double Beta(double* i, double* j) {
	double temp1 = 0.0, temp2 = 0.0;
	for(int k = 0; k < m; k++) {
		temp1 += i[k] * i[k];
		temp2 += j[k] * j[k];
	}
	return temp1 / temp2;
}

double** init(int n) {
	double** P = new double*[n+1];
	double* alpha = new double[n+1];
	double* beta = new double[n+1];
	beta[0] = 0.0;
	alpha[0] = 0.0;
	for(int i = 0; i <= n; i++)
		P[i] = new double[m];
	for(int i = 0; i < m; i++)
		P[0][i] = 1.0;
	for(int i = 1; i <= n; i++) {
		alpha[i] = Alpha(P[i - 1]);
		for(int j = 0; j < m; j++) {
			if(i==1)
				P[i][j] = (x[j] - alpha[i])*P[i - 1][j];
			else
				P[i][j] = (x[j] - alpha[i])*P[i - 1][j] - beta[i - 1] * P[i - 2][j];
		}
		beta[i] = Beta(P[i], P[i - 1]);
	}
	double *a = new double[n+1];
	for(int k = 0; k <= n; k++) {
		double temp1 = 0.0, temp2 = 0.0;
		for(int i = 0; i < m; i++) {
			temp1 += P[k][i] * y[i];
			temp2 += P[k][i] * P[k][i];
		}
		a[k] = temp1 / temp2;
	}
	double** arg=new double*[3];
	arg[0] = alpha;
	arg[1] = beta;
	arg[2] = a;
	return arg;
}

double S(double x, double** arg ,int n) {
	double P = 1.0, _P = 0.0,__P=0.0;
	double result = arg[2][0];
	for(int i = 1; i <= n; i++) {
		__P = _P;
		_P = P;
		P = (x - arg[0][i])*P - arg[1][i - 1] * __P;
		result += arg[2][i] * P;
	}
	return result;
}

int main() {
	ofstream file("data.txt");
	for(int n = 1; n <= 3; n++) {
		cout << n << endl;
		double** arg = init(n);
		for(double i = -2.0; i < 3.0; i += 0.02) {
			file << i << '\t' << S(i, arg, n) << endl;
		}
		for(int i = 0; i < m; i++) {
			cout << x[i] << '\t' << S(x[i], arg, n)-y[i] << endl;
		}
	}
	file.close();
	system("pause");
	return 0;
}