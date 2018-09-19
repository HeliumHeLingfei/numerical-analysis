#include <iostream>
#include <math.h>
using namespace std;

int main() {
	double ln2 = 0.693147190546, epsilon = 0.00005, result = 0, n = 1, m = -1, q = -1;
	while(abs(ln2 - result) > epsilon) {
		q = q * m;
		result += q / n;
		n++;
	}
	cout << "result:" << n << endl;
	system("pause");
	return 0;
}