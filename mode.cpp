#include <iostream>
#include <cmath>
using namespace std;

int main ()
{
	double x = 14.86, intPart, fractPart;
	int k;
	fractPart = modf(x, &intPart);
    k = intPart;
	cout << x << " = " << k << " + " << fractPart << endl;
	
	x = -31.201;
	fractPart = modf(x, &intPart);
    k = intPart;
	cout << x << " = " << k << " + " << fractPart << endl;

	return 0;
}