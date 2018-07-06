#include "src/geometric_cutting.h"

#include "gtest/gtest.h"

namespace {

using namespace R2P;

TEST(GeometricCutting, getShiftedUnitCube) {
  RectangularCuboid shifted_cube = getShiftedUnitCube(-1.0, -1.0, -1.0);
  UnitCube unit_cube;
  for (int v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(shifted_cube.vertex_m[v].x_m,
                     unit_cube.vertex_m[v].x_m - 1.0);
    EXPECT_DOUBLE_EQ(shifted_cube.vertex_m[v].y_m,
                     unit_cube.vertex_m[v].y_m - 1.0);
    EXPECT_DOUBLE_EQ(shifted_cube.vertex_m[v].z_m,
                     unit_cube.vertex_m[v].z_m - 1.0);
  }
}

TEST(GeometricCutting, getCutTet) {}

}  // namespace
