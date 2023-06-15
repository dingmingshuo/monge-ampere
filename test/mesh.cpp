#include <gtest/gtest.h>
#include "ma.hpp"

TEST(MeshTest, PointLocationTest) {
  auto mesh = Mesh("./mesh/2e-2.obj");
  int id;
  id = mesh.point_location(Point(0.4, 0.5));
  EXPECT_EQ(id, 32);
  id = mesh.point_location(Point(1, 0.6));
  EXPECT_EQ(id, 29);
  id = mesh.point_location(Point(1.2, 1.2));
  EXPECT_EQ(id, -1);
}
