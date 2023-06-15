#include "ma.hpp"

#include <cmath>
#include <cstdio>
#include <algorithm>

int main() {
    auto f = [](Point p) { return 2 * pow(2 -  p.norm2(), -2); };
    auto u = [](Point p) { return - sqrt(2 - p.norm2()); };                
    std::string filename = "./mesh/2e-8.obj";
    auto mesh = Mesh(filename);
    auto u_gt = MeshFunction(mesh);
    u_gt.init(u);
    real h = std::pow(2.0, -8.0);
    real delta = std::pow(h, 1.0 / 2.0);
    real theta = std::pow(h, 1.0 / 2.0);
    auto u_eps = two_scale::perron(f, u, mesh, delta, theta);
    printf("Error |u-u_eps|_inf:\n%e\n", (u_gt - u_eps).norm_inf());

    return 0;
}