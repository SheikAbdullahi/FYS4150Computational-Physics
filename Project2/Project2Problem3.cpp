#include <iostream>
#include <armadillo>

// Function to find the maximum off-diagonal element of a symmetric matrix
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

int main() {
    // Define the test matrix A
    arma::mat A = {{1, 0, 0, 0.5},
                   {0, 1, -0.7, 0},
                   {0, -0.7, 1, 0},
                   {0.5, 0, 0, 1}};

    // Call the function to find the maximum off-diagonal element
    int k, l;
    double max_val = max_offdiag_symmetric(A, k, l);
    
    // Print the result
    std::cout << "The maximum off-diagonal element is: " << max_val << " at indices (" << k << ", " << l << ").\n";
    
    return 0;
}
