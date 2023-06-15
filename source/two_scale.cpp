#include "ma.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace two_scale {
    std::vector<Point> generate_s_theta(real theta) {
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

    real T_epsilon(MeshFunction& u, int id, std::vector<Point> s_theta, real delta) {
        int D = (int)s_theta.size() / 4;

        real ret = INF;
        #pragma omp parallel for reduction(min:ret)
        for (int i = 0; i < (int)s_theta.size() - D; i++) {
            Point v0 = s_theta[i];
            Point v1 = s_theta[i + D];

            real diff0 = u.second_difference(id, v0, delta);

            real diff1 = u.second_difference(id, v1, delta);

            real t = (std::max(diff0, 0.0) * std::max(diff1, 0.0)) + 
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
        real rho;
    };

    std::vector<std::vector<PreparedSecondDiff>> prepare_for_second_difference(
        Mesh& mesh, std::vector<int>& N0, std::vector<Point>& s_theta, real delta) {
        std::vector<std::vector<PreparedSecondDiff>> ret;
        ret.resize(mesh.num_points);
        #pragma omp parallel for
        for (auto i: N0) {
            std::vector<PreparedSecondDiff> prepared;
            prepared.resize(s_theta.size());
            Point p = mesh.points[i];
            for (int j = 0; j < (int)s_theta.size(); j++) {
                Point v = s_theta[j];
                real rho = 1;

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

    real T_epsilon_prepared(MeshFunction& u, int id, std::vector<Point> s_theta, real delta,
                            std::vector<std::vector<PreparedSecondDiff>>& prepared) {
        int D = (int)s_theta.size() / 4;
        std::vector<real> diff;
        diff.resize((int)s_theta.size());
        for (int i = 0; i < (int)s_theta.size(); i++) {
            diff[i] = u.second_difference(id, prepared[id][i].p_plus, prepared[id][i].e_plus_id, 
                                             prepared[id][i].p_minus, prepared[id][i].e_minus_id, prepared[id][i].rho, delta);
        }

        real ret = INF;
        for (int i = 0; i < (int)s_theta.size() - D; i++) {
            real diff0 = diff[i];
            real diff1 = diff[i + D];
            real t = (std::max(diff0, 0.0) * std::max(diff1, 0.0)) + 
                     (std::min(diff0, 0.0) + std::min(diff1, 0.0));
            ret = std::min(ret, t);
        }

        return ret;
    }
#endif

    MeshFunction perron(real (*f)(Point), real (*g)(Point), Mesh& mesh, real delta, real theta, int p) {
        // Initialize super parameters
        real tolerance = 1e-5;
        real max_u = 100;
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
        real f_max = -INF;
        real g_min = INF;
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
        real pre_err = INF;
        for (int iter = 0; iter < max_iter; iter++) {
            for (auto i: N0) {
                real T_now = two_scale::T_eps(u, i, s_theta, delta, prepared);
                if (T_now > F[i]) {
                    real l = u[i];
                    real r = max_u;
                    while (std::fabs(r - l) > EPS) {
                        real mid = (l + r) / 2;
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
            real err = (T - F).norm2();
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

} // namespace perron
