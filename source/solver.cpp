#include "ma.hpp"

using namespace Eigen;

VectorXd solver::GMRES(const SparseMatrix<double> &A,
                       const VectorXd &b, int m, int max_iter,
                       double tolerance) {
    int n = (int)A.rows();
    VectorXd x = b;
    VectorXd r = b - A * x;
    MatrixXd H = MatrixXd::Zero(m + 1, m); // Hessenberg matrix
    MatrixXd V = MatrixXd::Zero(n, m);     // Othonormal vectors
    VectorXd w(n);

    double beta = r.norm();
    V.col(0) = r / beta;
    VectorXd e1_beta = VectorXd::Zero(m + 1);
    e1_beta(0) = beta;

    for (int iter = 0; iter < max_iter; iter++) {
        for (int j = 0; j < m; j++) {
            // Arnoldi process
            VectorXd Av = A * V.col(j);
            for (int i = 0; i <= j; i++) {
                H(i, j) = Av.dot(V.col(i));
            }

            // update V
            if (j != m - 1) {
                VectorXd v_hat = Av;
                for (int i = 0; i <= j; i++) {
                    v_hat -= H(i, j) * V.col(i);
                }
                H(j + 1, j) = v_hat.norm();
                V.col(j + 1) = v_hat / H(j + 1, j);
            }
        }

        // solve least squares problem and update x
        JacobiSVD<MatrixXd> svd(H, ComputeThinV | ComputeThinU);
        VectorXd y = svd.solve(e1_beta);
        x += V * y;

        // compute new residual
        r = b - A * x;
        beta = r.norm();
        if (beta < tolerance) {
            break;
        }
        V.col(0) = r / beta;
        e1_beta(0) = beta;
    }

    return x;
}
