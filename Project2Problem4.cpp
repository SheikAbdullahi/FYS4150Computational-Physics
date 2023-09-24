#include <iostream>
#include <armadillo>
#include <cmath>

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l) {
    double max_val = 0.0;
    for(int i = 0; i < A.n_rows; ++i) {
        for(int j = i + 1; j < A.n_cols; ++j) {
            double a_ij = fabs(A(i, j));
            if(a_ij > max_val) {
                max_val = a_ij;
                k = i;
                l = j;
            }
        }
    }
    return max_val;
}

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l) {
    double s, c;
    if(A(k, l) != 0.0) {
        double t, tau = (A(l, l) - A(k, k)) / (2.0 * A(k, l));
        if(tau >= 0)
            t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        else
            t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        c = 1 / sqrt(1 + t * t);
        s = c * t;
    } else {
        c = 1.0;
        s = 0.0;
    }
    
    double a_kk = A(k, k), a_ll = A(l, l);
    A(k, k) = c * c * a_kk - 2.0 * c * s * A(k, l) + s * s * a_ll;
    A(l, l) = s * s * a_kk + 2.0 * c * s * A(k, l) + c * c * a_ll;
    A(k, l) = 0.0;
    A(l, k) = 0.0;
    
    for(int i = 0; i < A.n_rows; ++i) {
        if(i != k && i != l) {
            double a_ik = A(i, k), a_il = A(i, l);
            A(i, k) = c * a_ik - s * a_il;
            A(k, i) = A(i, k);
            A(i, l) = c * a_il + s * a_ik;
            A(l, i) = A(i, l);
        }
        
        double r_ik = R(i, k), r_il = R(i, l);
        R(i, k) = c * r_ik - s * r_il;
        R(i, l) = c * r_il + s * r_ik;
    }
}

int main() {
    // Initialize the matrix A and the identity matrix R
    int N = 6;
    int n = N + 1;
    double h = 1.0 / n;
    
    // Define elements of matrix A
    double a = -1 / (h * h);
    double d = 2 / (h * h);
    
    // Initialize matrix A and matrix R
    arma::mat A(N, N, arma::fill::zeros);
    arma::mat R = arma::eye<arma::mat>(N, N);
    for(int i = 0; i < N; ++i) {
        A(i, i) = d;
        if(i < N - 1) {
            A(i, i + 1) = a;
            A(i + 1, i) = a;
        }
    }
    
    int iterations = 0;
    double epsilon = 1.0e-8;
    int max_iterations = N * N * N;
    double max_offdiag = max_offdiag_symmetric(A, k, l);
    while(fabs(max_offdiag) > epsilon && iterations < max_iterations) {
        int k, l;
        max_offdiag = max_offdiag_symmetric(A, k, l);
        jacobi_rotate(A, R, k, l);
        iterations++;
    }
    
    // Check for convergence
    if(iterations < max_iterations) {
        std::cout << "Jacobi method converged in " << iterations << " iterations.\n";
    } else {
        std::cerr << "Error: Jacobi method did not converge within " << max_iterations << " iterations.\n";
        return 1;
    }
    
    // Output the results
    arma::vec eigenvalues = A.diag();
    std::cout << "Eigenvalues: " << eigenvalues.t() << "\n";
    std::cout << "Eigenvectors: \n" << R << "\n";
    
    return 0;
}
