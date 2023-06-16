#include "ma.hpp"

#include <cmath>
#include <cstdio>
#include <algorithm>

int main() {
    auto f = [](Point p) { return std::max(1.0 - 0.2 / (p - Point(0.5, 0.5)).norm(), 0.0); };
    auto u = [](Point p) { return 0.5 * std::max((p - Point(0.5, 0.5)).norm() - 0.2, 0.0) *
                                        std::max((p - Point(0.5, 0.5)).norm() - 0.2, 0.0); };                
    std::string filename = "./mesh/2e-8.obj";
    auto mesh = Mesh(filename);
    auto u_gt = MeshFunction(mesh);
    u_gt.init(u);
    ma_real h = std::pow(2.0, -8.0);
    ma_real delta = std::pow(h, 4.0 / 5.0);
    ma_real theta = std::pow(h, 2.0 / 5.0);
    auto u_eps = two_scale::newton(f, u, mesh, delta, theta, 48);
    printf("Error |u-u_eps|_inf:\n%e\n", (u_gt - u_eps).norm_inf());

    return 0;
}