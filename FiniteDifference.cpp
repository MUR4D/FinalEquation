#include <iostream>
#include <vector>
#include <math.h>
#include "inv_matrix.h"
#include "progonka.h"


int main() {
	int n = 600, i;
	double alpha = 1, default_value = 0, h = 1.0 / (n - 1);
	std::vector <double> u(n, default_value);
	std::vector <double> u_n(n, default_value);
	std::vector <double> u_nk(n, default_value);

	u_nk = progonka(n, h, alpha);

	for (i = 0; i < n; i++)std::cout << "u[" << i << "] = " << u_nk[i] << std::endl;
	std::cout << std::endl;


	return 0;
}

