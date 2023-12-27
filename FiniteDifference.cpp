#include <iostream>
#include <vector>
#include <math.h>

std::vector<double> progonka(int n, double h, double alpha);
std::vector<std::vector<double>> inv_matr(int n, std::vector<std::vector<double>> matr);

int main() {
    int n = 11, i;
    double alpha = 1, default_value = 0, h = 1.0 / (n - 1);
    std::vector<double> u(n, default_value);
    std::vector<double> u_n(n, default_value);
    std::vector<double> u_nk(n, default_value);

    u_nk = progonka(n, h, alpha);

    for (i = 0; i < n; i++)
        std::cout << "u[" << i << "] = " << u_nk[i] << std::endl;
    std::cout << std::endl;

    return 0;
}

std::vector<double> progonka(int n, double h, double alpha) {
    int i, j, k1 = 0;
    double default_value = 0, sum = 0, eps = 1e-6, sum1;
    std::vector<std::vector<double>> A(n, std::vector<double>(n, default_value));
    std::vector<std::vector<double>> inv_A(n, std::vector<double>(n, default_value));
    std::vector<double> u_n(n, default_value);
    std::vector<double> u_nk(n, default_value);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j)
                A[i][j] = -2 / pow(h, 2);
            else if (i == j + 1)
                A[i][j] = 1 / pow(h, 2) - sinh((i + 1) * h) / (2 * h);
            else if (i == j - 1)
                A[i][j] = 1 / pow(h, 2) + sinh(i * h) / (2 * h);
            else
                A[i][j] = 0;
        }
    }

    A[0][0] = A[n - 1][n - 1] = 1;
    A[0][1] = A[n - 1][n - 2] = 0;

    inv_A = inv_matr(n, A);

    u_n[0] = -1;
    u_n[n - 1] = 1;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            u_nk[i] += inv_A[i][j] * u_n[j];
        }
    }

    while (k1 < 100) {
        for (i = 0; i < n; i++) {
            u_n[i] = -pow((i * h), 2) * u_nk[i] - pow(i * h, 2);
        }
        u_n[0] = 1;
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
        if (sum < eps)
            break;
        else {
            sum = 0;
        }
        k1++;
    }
    k1 = 0;

    std::cout << "sum = " << sum1 << std::endl << std::endl;

    return u_nk;
}

std::vector<std::vector<double>> inv_matr(int n, std::vector<std::vector<double>> matr) {
    int i, j, m;
    double del, max_a, b;
    std::vector<std::vector<double>> ematr(n, std::vector<double>(n, 0.0));

    for (i = 0; i < n; i++) {
        ematr[i][i] = 1.0;
    }

    for (m = 0; m < n; m++) {
        max_a = matr[m][m];

        for (i = m + 1; i < n; i++) {
            if (fabs(max_a) < fabs(matr[i][m])) {
                max_a = matr[i][m];
                std::swap(matr[m], matr[i]);
                std::swap(ematr[m], ematr[i]);
            }
        }

        for (i = m + 1; i < n; i++) {
            del = matr[i][m] / matr[m][m];
            for (j = 0; j < n; j++) {
                matr[i][j] -= del * matr[m][j];
                ematr[i][j] -= del * ematr[m][j];
            }
        }
    }

    for (m = n - 1; m > 0; m--) {
        for (i = m - 1; i >= 0; i--) {
            del = matr[i][m] / matr[m][m];
            for (j = 0; j < n; j++) {
                matr[i][j] -= del * matr[m][j];
                ematr[i][j] -= del * ematr[m][j];
            }
        }
    }

    for (i = 0; i < n; i++) {
        del = matr[i][i];
        for (j = 0; j < n; j++) {
            matr[i][j] /= del;
            ematr[i][j] /= del;
        }
    }

    return ematr;
}
