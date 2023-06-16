#include <gtest/gtest.h>
#include "ma.hpp"

TEST(SolverTest, GMRESTest) {
    Eigen::SparseMatrix<double> A(5, 5);

    // Fill in the sparse matrix with non-zero elements
    A.insert(0, 0) = 2.0;
    A.insert(0, 1) = -1.0;
    A.insert(1, 0) = -1.0;
    A.insert(1, 1) = 2.0;
    A.insert(1, 2) = -1.0;
    A.insert(2, 1) = -1.0;
    A.insert(2, 2) = 2.0;
    A.insert(2, 3) = -1.0;
    A.insert(3, 2) = -1.0;
    A.insert(3, 3) = 2.0;
    A.insert(3, 4) = -1.0;
    A.insert(4, 3) = -1.0;
    A.insert(4, 4) = 2.0;

    // Compress the sparse matrix
    A.makeCompressed();

    // Define the right-hand side vector
    VectorXd b(5);
    b << 1, 2, 3, 4, 5;

    // Initialize the solution vector with all zeros
    VectorXd x(5);
    x.setZero();

    // Define the tolerance for the solution
    double tolerance = 1e-10;
    int max_iter = 10;
    int m = 3;

    // Call GMRES routine
    x = solver::GMRES(A, b, m, max_iter, tolerance);

    // Check the solution
    EXPECT_NEAR((A * x - b).norm(), 0.0, 1e-6);
}