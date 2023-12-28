#include<iostream>
#include<cmath>
#include<vector>
#include "inv_matrix.h"

std::vector <double> progonka(int n, double h, double alpha) {
	int i, j, k1 = 0;
	double default_value = 0, sum = 0, eps = 1e-6, sum1;
	std::vector <std::vector<double>> A(n, std::vector<double>(n, default_value));
	std::vector <std::vector<double>> inv_A(n, std::vector<double>(n, default_value));
	std::vector <double> u_n(n, default_value);
	std::vector <double> u_nk(n, default_value);


	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)A[i][j] = -2.0 / pow(h, 2) + h * i;
			else if (i == j + 1)A[i][j] = 1.0 / pow(h, 2);
			else if (i == j - 1)A[i][j] = 1.0 / pow(h, 2);
			else A[i][j] = 0;
		}
	}

	for (j = 0; j < n; j++) {
		A[0][j] = 0;
		A[n - 1][j] = 0;
	}

	A[0][0] = -1 / (2 * h);
	A[0][2] = 1 / (2 * h);
	A[n - 1][n - 3] = -1.0 / (2 * h);
	A[n - 1][n - 1] = 1.0 / (2 * h);

	inv_A = inv_matr(n, A);

	for (i = 0; i < n; i++) {
		u_n[i] = sinh(h * i) - i * h;
	}
	u_n[0] = -1;
	u_n[n - 1] = 1;

	/*
	for(i = 0; i < n; i++){
		std::cout << u_n[i] << " ";
	}
	std::cout << std::endl;
	*/

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			u_nk[i] += inv_A[i][j] * u_n[j];
		}
	}

	/*
	for(i = 0; i < n; i++)std::cout << u_nk[i] << " ";
	std::cout << std::endl;
	*/

	while (k1 < 100) {
		for (i = 1; i < n - 1; i++) {
			u_n[i] = (((-sinh(i * h) * (u_nk[i + 1] - u_nk[i - 1])) / (2 * h)) - i * h);

			//u_n[i] = -sinh(i * h) * (u_nk[i + 1] - u_nk[i - 1])/(2*h)-i*h;
			//u_n[i] = -pow((i * h), 2) * u_nk[i] - pow(i * h, 2);
		}

		u_n[0] = -1;
		u_n[n - 1] = 1;

		for (i = 0; i < n; i++) {
			u_nk[i] = 0;
			for (j = 0; j < n; j++) {
				u_nk[i] += inv_A[i][j] * u_n[j];
			}
		}

		for (i = 0; i < n - 1; i++) {
			sum += h * pow((u_nk[i + 1] - u_nk[i]), 2);
			sum1 = sum;
		}
		//std::cout << "sum = " << sum << std::endl;
		if (sum < eps)break;
		else {
			sum = 0;
		}
		k1++;
	}
	k1 = 0;

	std::cout << "sum = " << sum1 << std::endl << std::endl;

	return u_nk;
}
