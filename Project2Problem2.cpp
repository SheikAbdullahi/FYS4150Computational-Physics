
#include <iostream>
#include <armadillo>

int main() {
    // Define the size of the matrix
    int N = 6;
    int n = N + 1;
    double h = 1.0 / n; // stepsize
    
    // Define elements of matrix A
    double a = -1 / (h * h);
    double d = 2 / (h * h);
    
    // Initialize matrix A
    arma::mat A(N, N, arma::fill::zeros);
    for(int i = 0; i < N; ++i) {
        A(i, i) = d;
        if(i < N - 1) {
            A(i, i + 1) = a;
            A(i + 1, i) = a;
        }
    }
    
    // Solve the eigenvalue problem using Armadillo
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    
    // Display the computed eigenvalues and compare with the analytical solution
    std::cout << "Computed Eigenvalues: " << eigval.t() << std::endl;
    for(int j = 1; j <= N; ++j) {
        double lambda_analytical = d + 2 * a * cos(j * M_PI / (N + 1));
        std::cout << "Analytical Eigenvalue for j=" << j << ": " << lambda_analytical << std::endl;
    }
    
    // Check eigenvectors (normalize and compare with analytical solution)
    for(int j = 1; j <= N; ++j) {
        arma::vec v_computed = eigvec.col(j - 1);
        v_computed = arma::normalise(v_computed); // Normalize the eigenvector
        
        // Print and compare with analytical eigenvectors
        std::cout << "Computed Eigenvector for j=" << j << ": " << v_computed.t() << std::endl;
        
        arma::vec v_analytical(N);
        for(int i = 1; i <= N; ++i) {
            v_analytical(i - 1) = sin(i * j * M_PI / (N + 1));
        }
        v_analytical = arma::normalise(v_analytical); // Normalize the analytical eigenvector
        std::cout << "Analytical Eigenvector for j=" << j << ": " << v_analytical.t() << std::endl;
    }
    
    return 0;
}
