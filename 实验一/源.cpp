#include <iostream>
#include <fstream>

using namespace std;

double* d;
double* M;

double f(double x) {
	double i = 1.0 + 16.0 * x*x;
	return 1.0 / i;
}

double l(double x, double k, double n) {
	double h = 10.0 / n;
	double xk = -5.0 + k * h;
	double w1 = 1.0, w2 = 1.0;
	for(double i = -5.0; i <= 5.0; i += h) {
		if(i == xk)
			continue;
		w1 = w1 * (x - i);
		w2 = w2 * (xk - i);
	}
	return w1 / w2;
}

double L(double x, double n) {
	double result = 0;
	double h = 10.0 / n;
	for(double k = 0.0; k <= n; k++) {
		result += f(-5 + k * h)*l(x, k, n);
	}
	return result;
}

double spl2(double x[], double y[], int n, double ddy1, double ddyn, double t[], int m, double z[]) {
	int i, j;
	double h0, h1, alpha, beta, *s, *dy;
	s = new double[n];
	dy = new double[n];
	dy[0] = -0.5;
	h0 = x[1] - x[0];
	s[0] = 3.0*(y[1] - y[0]) / (2.0*h0) - ddy1 * h0 / 4.0;
	for(j = 1; j <= n - 2; j++) {
		h1 = x[j + 1] - x[j];
		alpha = h0 / (h0 + h1);
		beta = (1.0 - alpha)*(y[j] - y[j - 1]) / h0;
		beta = 3.0*(beta + alpha * (y[j + 1] - y[j]) / h1);
		dy[j] = -alpha / (2.0 + (1.0 - alpha)*dy[j - 1]);
		s[j] = (beta - (1.0 - alpha)*s[j - 1]);
		s[j] = s[j] / (2.0 + (1.0 - alpha)*dy[j - 1]);
		h0 = h1;
	}
	dy[n - 1] = (3.0*(y[n - 1] - y[n - 2]) / h1 + ddyn * h1 / 2.0 - s[n - 2]) / (2.0 + dy[n - 2]);
	for(j = n - 2; j >= 0; j--)        dy[j] = dy[j] * dy[j + 1] + s[j];
	for(j = 0; j <= n - 2; j++)        s[j] = x[j + 1] - x[j];
	for(j = 0; j <= m - 1; j++) {
		if(t[j] >= x[n - 1]) i = n - 2;
		else {
			i = 0;
			while(t[j]>x[i + 1]) i = i + 1;
		}
		h1 = (x[i + 1] - t[j]) / s[i];
		h0 = h1 * h1;
		z[j] = (3.0*h0 - 2.0*h0*h1)*y[i];
		z[j] = z[j] + s[i] * (h0 - h0 * h1)*dy[i];
		h1 = (t[j] - x[i]) / s[i];
		h0 = h1 * h1;
		z[j] = z[j] + (3.0*h0 - 2.0*h0*h1)*y[i + 1];
		z[j] = z[j] - s[i] * (h0 - h0 * h1)*dy[i + 1];
	}
	delete[] s;
	delete[] dy;
	return 0.0;
}

int main() {
	double a = f(-5) - L(-5, 10);
	double x1[11];
	double x2[21];
	double y1[11];
	double y2[21];
	double t[10001];
	double z1[1001], z2[1001];
	ofstream file("data.txt");
	int k = 0;
	for(double i = -5.0; i <= 5.0; i += 0.01) {
		t[k++] = i;
		file << i << '\t' << f(i) << endl;
	}
	file << "L10________________________" << endl;
	for(double i = -5.0; i <= 5.0; i += 0.01) {
		file << i << '\t' << L(i, 10.0) << endl;
	}
	file << "L20************************" << endl;
	for(double i = -5.0; i <= 5.0; i += 0.01) {
		file << i << '\t' << L(i, 20.0) << endl;
	}
	for(int i = 0; i <= 10; i++) {
		x1[i] = -5.0 + (double)i * 1;
		y1[i] = f(x1[i]);
	}
	for(int i = 0; i <= 20; i++) {
		x2[i] = -5.0 + (double)i * 0.5;
		y2[i] = f(x2[i]);
	}
	double m = 0.0005950261379281692;
	spl2(x1, y1, 10, m, m, t, 1001, z1);
	spl2(x2, y2, 20, m, m, t, 1001, z2);
	file << "S10************************" << endl;
	for(int i = 0; i <= 1000; i++) {
		file << t[i] << '\t' << z1[i] << endl;
	}
	file << "S20************************" << endl;
	for(int i = 0; i <= 1000; i++) {
		file << t[i] << '\t' << z2[i] << endl;
	}

	cout<<"f-L10:" << f(4.8) - L(4.8, 10) << endl;
	cout << "f-L20:" << f(4.8) - L(4.8, 20) << endl;
	cout << "f-S10:" << f(4.8) - z1[980] << endl;
	cout << "f-S20:" << f(4.8) - z2[980] << endl;
	file.close();
	system("pause");
	return 0;
}