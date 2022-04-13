/* test/abscissae.cpp
 * Test of constructor-algorithms of GaussKronrod class.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "gauss-kronrod.hpp"

int
main (int argc, char** argv)
{
	using namespace std;
	size_t m = (argc > 1) ? atoi(argv[1]) : 10;	// Legendre degree
	
	GaussKronrod<long double> GK(m);
	// machine epsilon ~ 10^{-d} :
	double eps_dbl = (double) GK.get_eps();
	int digits = ceil(fabs(log10(eps_dbl)));
	
	cout << "epsilon: " << GK.get_eps() << endl;
	cout << "Abscissae:" << endl;
	for (int k = 0; k < m+1; ++k) {
		cout.width(digits+5); cout.precision(digits);
		cout << GK.xgk(k) << endl;
	}
	cout << "Gauss-Legendre Weights:" << endl;
	for (int k = 0; k < (m+1)/2; ++k) {
		cout.width(digits+5); cout.precision(digits);
		cout << GK.wg(k) << endl;
	}
	cout << "Gauss-Kronrod Weights:" << endl;
	for (int k = 0; k < m+1; ++k) {
		cout.width(digits+5); cout.precision(digits);
		cout << GK.wgk(k) << endl;
	}
	
	return 0;
}
