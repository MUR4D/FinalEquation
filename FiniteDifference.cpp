#include <iostream>
#include <cmath>
#include <vector>

std::vector<std::vector<double>> matrix_A(int n, double h) {
    std::vector<std::vector<double>> matrix;

    for (int i = 0; i != n; i++) {
        std::vector<double> v;
        for (int j = 0; j != n; j++) {
            v.push_back(0.0);
        }
        matrix.push_back(v);
    }

    for (int i = 0; i < n; i++) {
        double x_i = i * h;
        double sinh_xi = sinh(x_i);

        matrix[i][i] = -2.0 - (x_i * h * h + sinh_xi) / (h * h); 

        if (i > 0) {
            matrix[i][i - 1] = 1.0 - sinh_xi / (2.0 * h); 
        }
        if (i < n - 1) {
            matrix[i][i + 1] = 1.0 + sinh_xi / (2.0 * h);
        }
    }

    
    matrix[0][0] = 1.0 - sinh(0) / (2.0 * h); 
    matrix[0][1] = 1.0 + sinh(0) / (2.0 * h);
    matrix[n - 1][n - 2] = 1.0 - sinh(1) / (2.0 * h);
    matrix[n - 1][n - 1] = 1.0 + sinh(1) / (2.0 * h); 
        
    return matrix;
}

std::vector<double> solveGauss(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = b.size();

    for (int i = 0; i < n; i++) {
        double maxEl = fabs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > maxEl) {
                maxEl = fabs(A[k][i]);
                maxRow = k;
            }
        }

        for (int k = i; k < n; k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }
        double tmp = b[maxRow];
        b[maxRow] = b[i];
        b[i] = tmp;

        for (int k = i + 1; k < n; k++) {
            double coeff = -A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                if (i == j) {
                    A[k][j] = 0.0;
                }
                else {
                    A[k][j] += coeff * A[i][j];
                }
            }
            b[k] += coeff * b[i];
        }
    }

    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; k--) {
            b[k] -= A[k][i] * x[i];
        }
    }

    return x;
}

std::vector<double> vector_b(int n, double h) {
    std::vector<double> vec(n, 0.0);

    for (int i = 0; i < n; i++) {
        double x_i = i * h;
        vec[i] = -x_i * (x_i + 1) * h * h; 
    }

   
    vec[0] -= (1.0 - sinh(0) / (2.0 * h)); 
    vec[n - 1] -= (1.0 + sinh(1) / (2.0 * h));

    return vec;
}

double dist(int n, double h, std::vector<double>& u, std::vector<double>& u1)
{
    double sum = 0.0;
    for (int i = 0; i != n; i++)
    {
        sum += h * (u1[i] - u[i]) * (u1[i] - u[i]);
    }
    return sum;
}
int main() {
    int n = 11;
    double h = 0.1;
    double epsilon = 1e-6; // Заданная точность

    std::vector<std::vector<double>> matrix = matrix_A(n, h);
    std::vector<double> vec = vector_b(n, h);

    std::vector<double> solution = solveGauss(matrix, vec);

    std::vector<double> initial_solution(n, 0.0); 

  

    std::vector<double> previous_solution = initial_solution;

    while (dist(n, h, solution, previous_solution) >= epsilon) {
        previous_solution = solution;
        vec = vector_b(n, h);
        solution = solveGauss(matrix, vec);
    }


    std::cout << "\nSolution:" << "\n";
    for (int i = 0; i < n; i++) {
        std::cout << "x[" << i << "]= " << solution[i] << "\n";
    }

    return 0;
}
