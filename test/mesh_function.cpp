#include <gtest/gtest.h>
#include "ma.hpp"
#include <cmath>

TEST(MeshFunctionTest, InterpolateQuadraticTest) {
    auto mesh = Mesh("./mesh/2e-5.obj");
    auto u = MeshFunction(mesh);
    auto f = [](Point p) { return p.x * p.x + p.y * p.y; };
    u.init(f);
    
    Point p(0.3, 0.5);
    auto delta = u.interpolate(p) - f(p);
    EXPECT_NEAR(std::fabs(delta), 0, 1e-3);

    p = Point(0.6, 0.7);
    delta = u.interpolate(p) - f(p);
    EXPECT_NEAR(std::fabs(delta), 0, 1e-3);
}

TEST(MeshFunctionTest, InterpolateLinearTest) {
    auto mesh = Mesh("./mesh/2e-5.obj");
    auto u = MeshFunction(mesh);
    auto f = [](Point p) { return 2 * p.x + 3 * p.y; };
    u.init(f);
    
    Point p(0.3, 0.5);
    auto delta = u.interpolate(p) - f(p);
    EXPECT_NEAR(std::fabs(delta), 0, 1e-3);

    p = Point(0.6, 0.7);
    delta = u.interpolate(p) - f(p);
    EXPECT_NEAR(std::fabs(delta), 0, 1e-3);
}

TEST(MeshFunctionTest, SecondDiffenceQuadraticTest) {
    auto mesh = Mesh("./mesh/2e-5.obj");
    auto u = MeshFunction(mesh);
    auto f = [](Point p) { return p.x * p.x + p.y * p.y; };
    u.init(f);
    
    real delta = std::sqrt(1.0 / 32.0);
    Point v(1, 0);
    auto diff = u.second_difference(1210, Point(0, 1), delta);
    EXPECT_NEAR(std::fabs(diff), 2, 0.1);
    diff = u.second_difference(1231, Point(1, 0), delta);
    EXPECT_NEAR(std::fabs(diff), 2, 0.1);
}