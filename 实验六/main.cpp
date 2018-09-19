#include <iostream>
#include <math.h>
#include <iomanip>

using namespace::std;
const double e= 1.0;
const double a= 0.5;
const int n= 101;
const double h= 0.01;
const double theta = 0.000000000001;

double labeled_y(double x) {//�õ�׼ȷֵ
	return (1.0 - a) * ((1.0 - exp(-x / e)) / (1.0 - exp(-1.0 / e))) + a * x;
}

double* label() {//�õ�׼ȷ����
	double* re = new double[n];
	double x = 0.0;
	for(int i = 0; i < n; i++) {
		re[i] = labeled_y(x);
		x += h;
	}
	return re;
}

double** getA() {
	double** A = new double*[n];
	for(int i = 0; i < n; i++) {
		A[i] = new double[n];
	}
	double a1 = -(2 * e + h);
	double a2 = e + h;
	A[0][0] = a1, A[0][1] = a2;
	for(int i = 1; i < n - 1; i++) 
		A[i][i - 1] = e, A[i][i] = a1, A[i][i + 1] = a2;
	A[n - 1][n - 2] = e, A[n - 1][n - 1] = a1;
	return A;
}

double* getb() {
	double* b = new double[n];
	for(int i = 0; i < n; i++)
		b[i] = a * h*h;
	return b;
}

double** init_x() {//��ʼ��x
	double **x = new double*[2];
	x[0] = new double[n];
	x[1] = new double[n];
	for(int i = 0; i < n; i++) {
		x[0][i] = 0.0;
		x[1][i] = 1.0;
	}
	x[0][n - 1] = 1.0;//�߽�����
	return x;
}

bool restrain(double* x1, double* x2) {//�ж��Ƿ�������
	double tem = 0;
	for(int i = 0; i < n; i++)
		if(abs(x1[i] - x2[i]) > theta) 
			return false;
	return true;
}

double* Jacobi(double** a,double* b) {
	double **x = init_x();
	int turn = 0, nturn = 1;//������������������x(k),x(k+1)�����Ϳռ临�Ӷ�
	int num = 0;
	while(!restrain(x[0], x[1])) {
		num++;
		turn = !turn;
		nturn = !turn;
		x[turn][0] = 0;
		for(int i = 1; i < n-1; i++) 
			x[turn][i] = (b[i] - a[i][i - 1] * x[nturn][i - 1] - a[i][i + 1] * x[nturn][i + 1]) / a[i][i];
		x[turn][n-1]= 1;
	}
	cout << num << endl;
	return x[turn];
}

double* GC(double** a, double* b) {
	double **x = init_x();
	int turn = 0, nturn = 1;
	int num = 0;
	while(!restrain(x[0], x[1])) {
		num++;
		turn = !turn;
		nturn = !turn;
		x[turn][0] = 0;//���ֱ߽�����
		for(int i = 1; i < n - 1; i++)
			x[turn][i] = (b[i] - a[i][i - 1] * x[turn][i - 1] - a[i][i + 1] * x[nturn][i + 1]) / a[i][i];
		x[turn][n - 1] = 1;
	}
	cout << num << endl;
	return x[turn];
}

double* SOR(double** a, double* b, double w) {
	double **x = init_x();
	int turn = 0, nturn = 1;
	int num = 0;
	while(!restrain(x[0], x[1])) {
		num++;
		turn = !turn;
		nturn = !turn;
		x[turn][0] = 0;
		for(int i = 1; i < n - 1; i++)
			x[turn][i] = w * ((b[i] - a[i][i - 1] * x[turn][i - 1] - a[i][i + 1] * x[nturn][i + 1]) / a[i][i]) + (1 - w) * x[nturn][i];//��Ȩƽ��
		x[turn][n - 1] = 1;
	}
	cout << num << endl;
	return x[turn];
}

int main() {
	double **a = getA();
	double *b = getb();
	double* y = label();
	double* x1 = Jacobi(a,b);
	double* x2 = GC(a,b);
	SOR(a, b, 1.0);
	SOR(a, b, 1.2);
	SOR(a, b, 1.4);
	SOR(a, b, 1.6);
	SOR(a, b, 1.8);
	SOR(a, b, 1.9);
	double* x3 = SOR(a, b, 1.94);
	SOR(a, b, 1.95);
	double q = 0.0;
	for(int i = 0; i < n; i++) {
		//cout << setprecision(4) << y[i] << ' ' << x1[i] << " " << x2[i]<<" " << x3[i] << endl;//������
		q += 0.01;
		int num = 1;
		if(abs(x1[i] - y[i]) < abs(x2[i] - y[i])) {
			if(abs(x1[i] - y[i]) < abs(x3[i] - y[i]))
				num = 1;
			else
				num = 3;
		}
		else {
			if(abs(x2[i] - y[i]) < abs(x3[i] - y[i]))
				num = 2;
			else
				num = 3;
		}
		int max = 1;
		if(abs(x1[i] - y[i]) > abs(x2[i] - y[i])) {
			if(abs(x1[i] - y[i]) > abs(x3[i] - y[i]))
				max = 1;
			else
				max = 3;
		}
		else {
			if(abs(x2[i] - y[i]) > abs(x3[i] - y[i]))
				max = 2;
			else
				max = 3;
		}
		printf("%d,%d: ", num,max);//�ж���ӽ�����Զ����
		cout << setprecision(4) <<  x1[i]- y[i] << " " << x2[i]- y[i] << " " << x3[i]- y[i] << endl;//������
	}
	int n;
	cin >> n;
}