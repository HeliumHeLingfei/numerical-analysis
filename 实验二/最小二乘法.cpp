#include <iostream>
#include <fstream>

using namespace std;

int m = 7;
double x1[7] = { -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0 };
double y[7] = { -4.467,-0.452, 0.551, 0.048, -0.447, 0.549, 4.552 };

double* gaussin(double** a, double* b,int n) {
	double* c=new double[n];  
	for(int k = 0; k < n - 1; k++) {
		for(int i = k + 1; i < n; i++)
			c[i] = a[i][k] / a[k][k];
		for(int i = k + 1; i < n; i++) {
			for(int j = 0; j < n; j++) {
				a[i][j] = a[i][j] - c[i] * a[k][j];
			}
			b[i] = b[i] - c[i] * b[k];
		}
	}

	double* x=new double[n];
	x[n - 1] = b[n - 1] / a[n - 1][n - 1];
	for(int i = n - 2; i >= 0; i--) {
		double sum = 0;
		for(int j = i + 1; j < n; j++) {
			sum += a[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / a[i][i];
	}
	return x;
}

double* init(int n) {
	double** G = new double*[n+1];
	double** x = new double*[n+1];
	double* b = new double[n+1];
	for(int i = 0; i <= n; i++) {
		G[i] = new double[n+1];
		x[i] = new double[m];
	}
	for(int k = 0; k < m; k++) {
		x[0][k] = 1;
		x[1][k] = x1[k];
		for(int i = 2; i <= n; i++) 
			x[i][k] = x[i - 1][k] * x1[k];
	}
	for(int i = 0; i <= n; i++) {
		for(int j = 0; j <= n; j++) {
			double sum = 0.0;
			for(int k = 0; k < m; k++)
				sum += x[i][k] * x[j][k];
			G[i][j] = sum;
		}
		double sum = 0.0;
		for(int k = 0; k < m; k++)
			sum += x[i][k] *y[k];
		b[i] = sum;
	}
	double* a = gaussin(G, b, n+1);
	return a;
}

double S(double* a,double _x,int n) {
	double x = 1.0,result=0.0;
	for(int i = 0; i <= n; i++) {
		result += x * a[i];
		x *= _x;
	}
	return result;
}

int main() {
	ofstream file("data2.txt");
	for(int n = 1; n <= 3; n++) {
		double* a = init(n);
		cout << "y = ";
		for(int i = 0; i <= n; i++) {
			cout << a[i] << " * x^" << i;
			if(i != n)
				cout << " + ";
		}
		cout << endl;
		for(double x =-2.0 ; x < 3.0; x+=0.02) {
			file <<x<<'\t'<< S(a, x, n) << endl;
		}
		for(int i = 0; i < m; i++) {
			cout << x1[i] << ": " << S(a, x1[i], n)-y[i] << endl;
		}
	}
	
	int k;
	cin >> k;
	return 0;
}