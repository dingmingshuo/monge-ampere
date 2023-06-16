#include "ma.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <map>
#include <iostream>

namespace two_scale {
    std::vector<Point> generate_s_theta(ma_real theta) {
        // PI / (2.0 * std::acos(theta / 2.0) is a calculated minimum of p.
        // ceil(x / 4) * 4 makes it possible to have orthogonal pairs.
        int p = (int)std::ceil((PI / (std::acos(1 - theta * theta / 2.0))) / 4.0) * 4;
        std::vector<Point> s_theta;

        for (int k = 0; k < p; k++) {
            s_theta.push_back(Point(std::cos(2.0 * k * PI / p),
                                    std::sin(2.0 * k * PI / p)));
        }

        return s_theta;
    }

    std::vector<Point> generate_s_theta_by_p(int p) {
        // PI / (2.0 * std::acos(theta / 2.0) is a calculated minimum of p.
        // ceil(x / 4) * 4 makes it possible to have orthogonal pairs.
        std::vector<Point> s_theta;

        for (int k = 0; k < p; k++) {
            s_theta.push_back(Point(std::cos(2.0 * k * PI / p),
                                    std::sin(2.0 * k * PI / p)));
        }

        return s_theta;
    }

    ma_real T_epsilon(MeshFunction& u, int id, std::vector<Point> s_theta, ma_real delta) {
        int D = (int)s_theta.size() / 4;

        ma_real ret = INF;
        #pragma omp parallel for reduction(min:ret)
        for (int i = 0; i < (int)s_theta.size() - D; i++) {
            Point v0 = s_theta[i];
            Point v1 = s_theta[i + D];

            ma_real diff0 = u.second_difference(id, v0, delta);

            ma_real diff1 = u.second_difference(id, v1, delta);

            ma_real t = (std::max(diff0, 0.0) * std::max(diff1, 0.0)) + 
                     (std::min(diff0, 0.0) + std::min(diff1, 0.0));
            ret = std::min(ret, t);
        }

        return ret;
    }

#if TWO_SCALE_PREPARE
    struct PreparedSecondDiff {
        Point p_plus;
        Point p_minus;
        int e_plus_id;
        int e_minus_id;
        ma_real rho;
    };

    std::vector<std::vector<PreparedSecondDiff>> prepare_for_second_difference(
        Mesh& mesh, std::vector<int>& N0, std::vector<Point>& s_theta, ma_real delta) {
        std::vector<std::vector<PreparedSecondDiff>> ret;
        ret.resize(mesh.num_points);
        #pragma omp parallel for
        for (auto i: N0) {
            std::vector<PreparedSecondDiff> prepared;
            prepared.resize(s_theta.size());
            Point p = mesh.points[i];
            for (int j = 0; j < (int)s_theta.size(); j++) {
                Point v = s_theta[j];
                ma_real rho = 1;

                if (v.x != 0) {
                    rho = std::min(rho, sgn(v.x) * (1 - p.x) / (delta * v.x));
                    rho = std::min(rho, sgn(v.x) * p.x / (delta * v.x));
                }
                if (v.y != 0) {
                    rho = std::min(rho, sgn(v.y) * (1 - p.y) / (delta * v.y));
                    rho = std::min(rho, sgn(v.y) * p.y / (delta * v.y));
                }

                Point p_plus = p + v * delta * rho;
                Point p_minus = p - v * delta * rho;

                int e_plus_id = mesh.point_location(p_plus);
                int e_minus_id = mesh.point_location(p_minus);

                if (e_plus_id == -1 || e_minus_id == -1) {
                    throw std::runtime_error("Point not in mesh! (Invalid rho or delta)");
                }

                prepared[j] = PreparedSecondDiff({p_plus, p_minus, e_plus_id, e_minus_id, rho});
            }
            ret[i].assign(prepared.begin(), prepared.end());
        }
        return ret;
    }

    ma_real T_epsilon_prepared(MeshFunction& u, int id, std::vector<Point> s_theta, ma_real delta,
                            std::vector<std::vector<PreparedSecondDiff>>& prepared) {
        int D = (int)s_theta.size() / 4;
        std::vector<ma_real> diff;
        diff.resize((int)s_theta.size());
        for (int i = 0; i < (int)s_theta.size(); i++) {
            diff[i] = u.second_difference(id, prepared[id][i].p_plus, prepared[id][i].e_plus_id, 
                                             prepared[id][i].p_minus, prepared[id][i].e_minus_id, prepared[id][i].rho, delta);
        }

        ma_real ret = INF;
        for (int i = 0; i < (int)s_theta.size() - D; i++) {
            ma_real diff0 = diff[i];
            ma_real diff1 = diff[i + D];
            ma_real t = (std::max(diff0, 0.0) * std::max(diff1, 0.0)) + 
                     (std::min(diff0, 0.0) + std::min(diff1, 0.0));
            ret = std::min(ret, t);
        }

        return ret;
    }
#endif

    MeshFunction perron(ma_real (*f)(Point), ma_real (*g)(Point), Mesh& mesh, ma_real delta, ma_real theta, int p) {
        // Initialize super parameters
        ma_real tolerance = 1e-5;
        ma_real max_u = 100;
        int max_iter = 5000;
        std::printf("delta: %lf, theta: %lf\n", delta, theta);

        // Initialize mesh, meshfunctions and s_theta
        auto N0 = mesh.N0();
        auto Nb = mesh.Nb();
        std::printf("|N0|: %lu, |Nb|: %lu\n", N0.size(), Nb.size());
        std::vector<Point> s_theta;
        if (p == -1) {
            s_theta = two_scale::generate_s_theta(theta);
        } else {
            std::printf("p: %d, generate s_theta by p.\n", p);
            s_theta = two_scale::generate_s_theta_by_p(p);
        }
        std::printf("|s_theta|: %lu\n", s_theta.size());
        MeshFunction u(mesh);
        MeshFunction T(mesh);
        MeshFunction F(mesh);
        F.init(f);
        for (auto i: Nb) {
            F[i] = 0;
            T[i] = 0;
        }

        
#if TWO_SCALE_PREPARE
        // Set OpenMP parallel
        omp_set_num_threads((int)s_theta.size());
        printf("Number of OpenMP threads: %d\n", omp_get_max_threads());
        // Prepare for second difference
        std::printf("Preparing for second difference...\n");
        omp_set_num_threads((int)s_theta.size());
        auto prepared_start = omp_get_wtime();
        auto prepared = two_scale::prepare_for_second_difference(mesh, N0, s_theta, delta);
        auto prepared_end = omp_get_wtime();
        std::printf("Preparation done. Time: %lf s\n", prepared_end - prepared_start);
#else
        // Set OpenMP parallel
        omp_set_num_threads((int)s_theta.size());
        printf("Number of OpenMP threads: %d\n", omp_get_max_threads());
#endif

        // Set initial value of u
        ma_real f_max = -INF;
        ma_real g_min = INF;
        for (auto i: N0) {
            f_max = std::max(f_max, F[i]);
        }
        auto q = [=](Point p){ return 0.5 * std::sqrt(f_max) * (p.x * p.x + p.y * p.y); };
        for (auto i: Nb) {
            u[i] = g(mesh.points[i]);
            g_min = std::min(g_min, u[i] - q(mesh.points[i]));
        }
        for (auto i: N0) {
            u[i] = g_min + q(mesh.points[i]);
        }


#if TWO_SCALE_PREPARE
    #define T_eps T_epsilon_prepared
#else
    #define T_eps T_epsilon
#endif

        // Iteration
        std::printf("Start iteration...\n");
        auto start = omp_get_wtime();
        ma_real pre_err = INF;
        for (int iter = 0; iter < max_iter; iter++) {
            for (auto i: N0) {
                ma_real T_now = two_scale::T_eps(u, i, s_theta, delta, prepared);
                if (T_now > F[i]) {
                    ma_real l = u[i];
                    ma_real r = max_u;
                    while (std::fabs(r - l) > EPS) {
                        ma_real mid = (l + r) / 2;
                        u[i] = mid;
                        T_now = two_scale::T_eps(u, i, s_theta, delta, prepared);
                        if (T_now > F[i]) {
                            l = mid;
                        } else {
                            r = mid;
                        }
                    }
                }
            }

            // Update T
            for (auto i: N0) {
                T[i] = two_scale::T_eps(u, i, s_theta, delta, prepared);
            }

            // Check convergence
            ma_real err = (T - F).norm2();
            std::printf("[iter %4d] truncation err = %.10f (%e)\n", iter, err, err);
            if (std::fabs(pre_err - err) < tolerance) {
                break;
            }
            pre_err = err;
        }
        auto end = omp_get_wtime();
        std::printf("Iteration done. Time: %lf s\n", end - start);

        return u;
    }
    
    struct PreparedSecondDiffForNewton {
        Point p_plus;
        Point p_minus;
        int e_plus_id;
        int x_plus_id[3];
        ma_real plus_lambda[3];
        int e_minus_id;
        int x_minus_id[3];
        ma_real minus_lambda[3];
        ma_real rho;
        ma_real delta_hat;
    };

    std::vector<std::vector<PreparedSecondDiffForNewton>> prepare_for_second_difference_for_newton(
        Mesh& mesh, std::vector<int>& N0, std::vector<Point>& s_theta, ma_real delta) {
        std::vector<std::vector<PreparedSecondDiffForNewton>> ret;
        ret.resize(mesh.num_points);
        #pragma omp parallel for
        for (auto i: N0) {
            std::vector<PreparedSecondDiffForNewton> prepared;
            prepared.resize(s_theta.size());
            Point p = mesh.points[i];
            for (int j = 0; j < (int)s_theta.size(); j++) {
                Point v = s_theta[j];
                ma_real rho = 1;

                if (v.x != 0) {
                    rho = std::min(rho, sgn(v.x) * (1 - p.x) / (delta * v.x));
                    rho = std::min(rho, sgn(v.x) * p.x / (delta * v.x));
                }
                if (v.y != 0) {
                    rho = std::min(rho, sgn(v.y) * (1 - p.y) / (delta * v.y));
                    rho = std::min(rho, sgn(v.y) * p.y / (delta * v.y));
                }

                Point p_plus = p + v * delta * rho;
                Point p_minus = p - v * delta * rho;

                int e_plus_id = mesh.point_location(p_plus);
                int e_minus_id = mesh.point_location(p_minus);

                if (e_plus_id == -1 || e_minus_id == -1) {
                    throw std::runtime_error("Point not in mesh! (Invalid rho or delta)");
                }

                prepared[j] = PreparedSecondDiffForNewton();
                prepared[j].p_plus = p_plus;
                prepared[j].p_minus = p_minus;
                prepared[j].e_plus_id = e_plus_id;
                prepared[j].e_minus_id = e_minus_id;
                for (int k = 0; k < 3; k++) {
                    prepared[j].x_plus_id[k] = mesh.elements[e_plus_id].v[k];
                    prepared[j].x_minus_id[k] = mesh.elements[e_minus_id].v[k];
                }
                mesh.get_barycentric_coordinated(e_plus_id, p_plus, prepared[j].plus_lambda);
                mesh.get_barycentric_coordinated(e_minus_id, p_minus, prepared[j].minus_lambda);
                prepared[j].rho = rho;
                prepared[j].delta_hat = delta * rho;
            }
            ret[i].assign(prepared.begin(), prepared.end());
        }
        return ret;
    }

    std::pair<ma_real, int> T_epsilon_prepared_for_newton(MeshFunction& u, int id, std::vector<Point> s_theta, ma_real delta,
                        std::vector<std::vector<PreparedSecondDiffForNewton>>& prepared) {
        int D = (int)s_theta.size() / 4;
        std::vector<ma_real> diff;
        diff.resize((int)s_theta.size());
        for (int i = 0; i < (int)s_theta.size(); i++) {
            diff[i] = u.second_difference(id, prepared[id][i].p_plus, prepared[id][i].e_plus_id, 
                                             prepared[id][i].p_minus, prepared[id][i].e_minus_id, prepared[id][i].rho, delta);
        }

        ma_real ret = INF;
        int k = 0;
        for (int i = 0; i < (int)s_theta.size() - D; i++) {
            ma_real diff0 = diff[i];
            ma_real diff1 = diff[i + D];
            ma_real t = (std::max(diff0, 0.0) * std::max(diff1, 0.0)) + 
                     (std::min(diff0, 0.0) + std::min(diff1, 0.0));
            if (t < ret) {
                ret = t;
                k = i;
            }
        }

        return std::make_pair(ret, k);
    }

    MeshFunction poisson_on_coarse_mesh(ma_real (*f)(Point), ma_real (*g)(Point), Mesh& mesh) {
        // Solve poisson equation on a uniform mesh with diagram equals to h

        // Create a uniform mesh
        Mesh coarse_mesh(std::pow(2.0, -5.0));
        MeshFunction u0(coarse_mesh);
        MeshFunction u(mesh);
        ma_real h = std::pow(2.0, -5.0) / std::sqrt(2.0);

        // Create function
        auto q = [=](Point p) { return std::sqrt(2 * f(p)); };

        // Generate matrix A
        int n = (int)std::ceil(1 / h);
        h = 1.0 / n;
        int N = (n - 1) * (n - 1);
        Eigen::SparseMatrix<ma_real> A(N, N);
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - 1; j++) {
                int id = i * (n - 1) + j;
                A.insert(id, id) = -4 / (h * h);
                if (i != 0) {
                    A.insert(id, id - (n - 1)) = 1 / (h * h);
                }
                if (i != n - 2) {
                    A.insert(id, id + (n - 1)) = 1 / (h * h);
                }
                if (j != 0) {
                    A.insert(id, id - 1) = 1 / (h * h);
                }
                if (j != n - 2) {
                    A.insert(id, id + 1) = 1 / (h * h);
                }
            }
        }
        A.makeCompressed();
        
        // Generate vector b
        Eigen::VectorXd b(N);
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - 1; j++) {
                int id = i * (n - 1) + j;
                b(id) = q(Point((i + 1) * h, (j + 1) * h));
                if (i == 0) {
                    b(id) -= g(Point(0, (j + 1) * h)) / (h * h);
                }
                if (i == n - 2) {
                    b(id) -= g(Point(1, (j + 1) * h)) / (h * h);
                }
                if (j == 0) {
                    b(id) -= g(Point((i + 1) * h, 0)) / (h * h);
                }
                if (j == n - 2) {
                    b(id) -= g(Point((i + 1) * h, 1)) / (h * h);
                }
            }
        }

        // Solve Au = b
        auto x = solver::CG(A, b, 100, 1e-6);

        // Assign x to u0
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - 1; j++) {
                u0[(i + 1) * (n + 1) + (j + 1)] = x(i * (n - 1) + j);
            }
        }
        for (auto i: coarse_mesh.Nb()) {
            u0[i] = g(coarse_mesh.points[i]);
        }

        // Interpolate u0 to u
        for (auto i: mesh.N0()) {
            u[i] = u0.interpolate(mesh.points[i]);
        }
        for (auto i: mesh.Nb()) {
            u[i] = g(mesh.points[i]);
        }

        return u;
    }

    MeshFunction newton(ma_real (*f)(Point), ma_real (*g)(Point), Mesh& mesh, ma_real delta, ma_real theta, int p) {
        // Initialize super parameters
        ma_real tolerance = 1e-8;
        int max_iter = 100;
        std::printf("delta: %lf, theta: %lf\n", delta, theta);

        // Initialize mesh, meshfunctions and s_theta
        auto N0 = mesh.N0();
        auto Nb = mesh.Nb();
        std::printf("|N0|: %lu, |Nb|: %lu\n", N0.size(), Nb.size());
        std::vector<Point> s_theta;
        if (p == -1) {
            s_theta = two_scale::generate_s_theta(theta);
        } else {
            std::printf("p: %d, generate s_theta by p.\n", p);
            s_theta = two_scale::generate_s_theta_by_p(p);
        }
        std::printf("|s_theta|: %lu\n", s_theta.size());
        MeshFunction T(mesh);
        MeshFunction F(mesh);
        F.init(f);
        for (auto i: Nb) {
            F[i] = 0;
            T[i] = 0;
        }

        // Set OpenMP parallel
        printf("Number of OpenMP threads: %d\n", omp_get_max_threads());
        // Prepare for second difference
        std::printf("Preparing for second difference...\n");
        omp_set_num_threads((int)s_theta.size());
        auto prepared_start = omp_get_wtime();
        auto prepared = two_scale::prepare_for_second_difference_for_newton(mesh, N0, s_theta, delta);
        auto prepared_end = omp_get_wtime();
        std::printf("Preparation done. Time: %lf s\n", prepared_end - prepared_start);

        // Set initial value of u
        printf("Solving Poisson equation on coarse mesh...\n");
        auto u = poisson_on_coarse_mesh(f, g, mesh);
        printf("Poisson equation solved.\n");

        // Generate map for N0
        std::map<int, int> N0_map;
        for (int i = 0; i < (int)N0.size(); i++) {
            N0_map[N0[i]] = i;
        }

        // Newton iteration
        std::printf("Newton iteration...\n");
        auto start = omp_get_wtime();
        ma_real pre_err = INF;
        for (int iter = 0; iter < max_iter; iter++) {
            // Compute Jacobian Matrix DT
            Eigen::SparseMatrix<double> DT(N0.size(), N0.size());
            std::map<int, double> DT_triplets[mesh.num_points];
            #pragma omp parallel for
            for (auto i: N0) {
                auto [_, mn] = T_epsilon_prepared_for_newton(u, i, s_theta, delta, prepared);
                auto prepared_1 = prepared[i][mn];
                auto prepared_2 = prepared[i][mn + (int)s_theta.size() / 4];
                auto second_diff_1 = u.second_difference(i, prepared_1.p_plus, prepared_1.e_plus_id, 
                                                         prepared_1.p_minus, prepared_1.e_minus_id,
                                                         prepared_1.rho, delta);
                auto second_diff_2 = u.second_difference(i, prepared_2.p_plus, prepared_2.e_plus_id, 
                                                         prepared_2.p_minus, prepared_2.e_minus_id,
                                                         prepared_2.rho, delta);
                ma_real h_1_plus = second_diff_1 > 0;
                ma_real h_1_minus = 1 - h_1_plus;
                ma_real h_2_plus = second_diff_2 > 0;
                ma_real h_2_minus = 1 - h_2_plus;
                second_diff_1 = std::max(second_diff_1, 0.0);
                second_diff_2 = std::max(second_diff_2, 0.0);
                DT_triplets[i][i] =
                    -2.0 / (prepared_1.delta_hat * prepared_1.delta_hat) * (h_1_plus * second_diff_2 + h_1_minus)
                    -2.0 / (prepared_2.delta_hat * prepared_2.delta_hat) * (h_2_plus * second_diff_1 + h_2_minus);
                for (int k = 0; k < 3; k++) {
                    // j == 1
                    DT_triplets[i][prepared_1.x_plus_id[k]] +=
                        prepared_1.plus_lambda[k] / (prepared_1.delta_hat * prepared_1.delta_hat) * (h_1_plus * second_diff_2 + h_1_minus);
                    DT_triplets[i][prepared_1.x_minus_id[k]] +=
                        prepared_1.minus_lambda[k] / (prepared_1.delta_hat * prepared_1.delta_hat) * (h_1_plus * second_diff_2 + h_1_minus);
                    // j == 2
                    DT_triplets[i][prepared_2.x_plus_id[k]] +=
                        prepared_2.plus_lambda[k] / (prepared_2.delta_hat * prepared_2.delta_hat) * (h_2_plus * second_diff_1 + h_2_minus);
                    DT_triplets[i][prepared_2.x_minus_id[k]] +=
                        prepared_2.minus_lambda[k] / (prepared_2.delta_hat * prepared_2.delta_hat) * (h_2_plus * second_diff_1 + h_2_minus);
                }
            }
            for (auto i: N0) {
                for (auto j: DT_triplets[i]) {
                    if (N0_map.find(j.first) != N0_map.end()) {
                        DT.insert(N0_map[i], N0_map[j.first]) = j.second;
                    }
                }
            }
            DT.makeCompressed();

            // Prepare for right hand side
            Eigen::VectorXd rhs(N0.size());
            #pragma omp parallel for
            for (auto i: N0) {
                auto [T_eps, _] = T_epsilon_prepared_for_newton(u, i, s_theta, delta, prepared);
                rhs(N0_map[i]) = F[i] - T_eps;
            }

            // Solve linear system
            auto w = solver::GMRES(DT, rhs, 100, 100, 1e-6);

            // Calculate damping parameter
            ma_real tau = 1;
            while (tau > tolerance) {
                MeshFunction u_new(mesh);
                for (auto i: Nb) {
                    u_new[i] = g(mesh.points[i]);
                }
                for (auto i: N0) {
                    u_new[i] = u[i] + tau * w(N0_map[i]);
                }
                for (auto i: N0) {
                    auto [T_eps, _] = T_epsilon_prepared_for_newton(u_new, i, s_theta, delta, prepared);
                    T[i] = T_eps;    
                }
                ma_real err = (T - F).norm2();
                if (err < pre_err) {
                    break;
                }
                tau /= 2;
            }

            // Update u
            for (auto i: N0) {
                u[i] += tau * w(N0_map[i]);
            }

            // Update T
            for (auto i: N0) {
                auto [T_eps, _] = T_epsilon_prepared_for_newton(u, i, s_theta, delta, prepared);
                T[i] = T_eps;    
            }

            // Check convergence
            ma_real err = (T - F).norm2();
            std::printf("[iter %4d] truncation err = %.10f (%e), tau = %.10f\n", iter, err, err, tau);
            if (std::fabs(pre_err - err) < tolerance || err < tolerance) {
                break;
            }
            pre_err = err;
        }
        auto end = omp_get_wtime();

        // Output
        std::printf("Newton iteration done. Time: %lf s\n", end - start);

        return u;
    }

} // namespace perron
