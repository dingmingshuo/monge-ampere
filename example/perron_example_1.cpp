#include "ma.hpp"

#include <cmath>
#include <cstdio>
#include <algorithm>

int main() {
    auto f = [](Point p) { return (1 + p.norm2()) * std::exp(p.norm2()); };
    auto u = [](Point p) { return std::exp(p.norm2() / 2); };
    std::string filename = "./mesh/2e-8.obj";
    auto mesh = Mesh(filename);
    auto u_gt = MeshFunction(mesh);
    u_gt.init(u);
    ma_real h = std::pow(2.0, -8.0);
    ma_real delta = std::pow(h, 1.0 / 2.0);
    ma_real theta = std::pow(h, 1.0 / 2.0);
    auto u_eps = two_scale::perron(f, u, mesh, delta, theta);
    printf("Error |u-u_eps|_inf:\n%e\n", (u_gt - u_eps).norm_inf());

    return 0;
}