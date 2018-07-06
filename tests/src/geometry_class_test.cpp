#include "gtest/gtest.h"
#include "src/geometry_classes.h"

namespace {

using namespace R2P;

TEST(GeometryClasses, Pt) {
  Pt pt0;
  pt0.setLoc(2.0, 1.0, 0.0);
  Pt pt1(0.0, 1.0, 2.0);

  Pt diff_pt = pt0 - pt1;
  EXPECT_DOUBLE_EQ(diff_pt.x_m, 2.0);
  EXPECT_DOUBLE_EQ(diff_pt.y_m, 0.0);
  EXPECT_DOUBLE_EQ(diff_pt.z_m, -2.0);

  Pt plus_pt = pt0 + pt1;
  EXPECT_DOUBLE_EQ(plus_pt.x_m, 2.0);
  EXPECT_DOUBLE_EQ(plus_pt.y_m, 2.0);
  EXPECT_DOUBLE_EQ(plus_pt.z_m, 2.0);

  Pt right_mult_pt = pt0 * 4.0;
  EXPECT_DOUBLE_EQ(right_mult_pt.x_m, 8.0);
  EXPECT_DOUBLE_EQ(right_mult_pt.y_m, 4.0);
  EXPECT_DOUBLE_EQ(right_mult_pt.z_m, 0.0);

  Pt left_mult_pt = 4.0 * pt0;
  EXPECT_DOUBLE_EQ(right_mult_pt.x_m, left_mult_pt.x_m);
  EXPECT_DOUBLE_EQ(right_mult_pt.y_m, left_mult_pt.y_m);
  EXPECT_DOUBLE_EQ(right_mult_pt.z_m, left_mult_pt.z_m);

  Pt divide_pt = pt0 / 2.0;
  EXPECT_DOUBLE_EQ(divide_pt.x_m, 1.0);
  EXPECT_DOUBLE_EQ(divide_pt.y_m, 0.5);
  EXPECT_DOUBLE_EQ(divide_pt.z_m, 0.0);
}

TEST(GeometryClasses, RectangularCubic) {
  RectangularCuboid rect_cubic;
  RectangularCuboid two_point_set_rect_cubic(Pt(-0.5, -0.5, -0.5),
                                             Pt(0.5, 0.5, 0.5));
  for (int v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].x_m,
                     two_point_set_rect_cubic.vertex_m[v].x_m);
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].y_m,
                     two_point_set_rect_cubic.vertex_m[v].y_m);
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].z_m,
                     two_point_set_rect_cubic.vertex_m[v].z_m);
  }

  rect_cubic.shiftInX(-1.0);
  rect_cubic.shiftInY(-1.0);
  rect_cubic.shiftInZ(-1.0);
  for (int v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].x_m,
                     two_point_set_rect_cubic.vertex_m[v].x_m - 1.0);
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].y_m,
                     two_point_set_rect_cubic.vertex_m[v].y_m - 1.0);
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].z_m,
                     two_point_set_rect_cubic.vertex_m[v].z_m - 1.0);
  }

  rect_cubic.shift(1.0, 1.0, 1.0);
  for (int v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].x_m,
                     two_point_set_rect_cubic.vertex_m[v].x_m);
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].y_m,
                     two_point_set_rect_cubic.vertex_m[v].y_m);
    EXPECT_DOUBLE_EQ(rect_cubic.vertex_m[v].z_m,
                     two_point_set_rect_cubic.vertex_m[v].z_m);
  }

  EXPECT_DOUBLE_EQ(rect_cubic.volume(), 1.0);
  Pt centroid = rect_cubic.centroid();
  EXPECT_DOUBLE_EQ(centroid.x_m, 0.0);
  EXPECT_DOUBLE_EQ(centroid.y_m, 0.0);
  EXPECT_DOUBLE_EQ(centroid.z_m, 0.0);

  UnitCube unit_cube;
  for (int v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(unit_cube.vertex_m[v].x_m,
                     two_point_set_rect_cubic.vertex_m[v].x_m);
    EXPECT_DOUBLE_EQ(unit_cube.vertex_m[v].y_m,
                     two_point_set_rect_cubic.vertex_m[v].y_m);
    EXPECT_DOUBLE_EQ(unit_cube.vertex_m[v].z_m,
                     two_point_set_rect_cubic.vertex_m[v].z_m);
  }
}

TEST(GeometryClasses, Plane) {
  // Currently this test is empty
}

TEST(GeometryClasses, GeometricMoments) {
  GeometricMoments init_moment;
  EXPECT_DOUBLE_EQ(init_moment.volume_m, 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid_m.x_m, 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid_m.y_m, 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid_m.z_m, 0.0);

  GeometricMoments nonzero_moment(0.25, Pt(-0.5, 1.0, -8.0));
  EXPECT_DOUBLE_EQ(nonzero_moment.volume_m, 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moment.centroid_m.x_m, -0.5);
  EXPECT_DOUBLE_EQ(nonzero_moment.centroid_m.y_m, 1.0);
  EXPECT_DOUBLE_EQ(nonzero_moment.centroid_m.z_m, -8.0);

  init_moment += nonzero_moment;
  EXPECT_DOUBLE_EQ(init_moment.volume_m, 0.25);
  EXPECT_DOUBLE_EQ(init_moment.centroid_m.x_m, -0.5);
  EXPECT_DOUBLE_EQ(init_moment.centroid_m.y_m, 1.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid_m.z_m, -8.0);
}

TEST(GeometryClasses, PhaseMoments) {
  GeometricMoments liquid_moments(0.25, Pt(-0.5, 1.0, -8.0));
  GeometricMoments gas_moments(0.80, Pt(0.5, -1.0, 8.0));

  PhaseMoments blank_moments;
  PhaseMoments nonzero_moments(liquid_moments, gas_moments);
  EXPECT_DOUBLE_EQ(nonzero_moments.liquid_m.volume_m, 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moments.liquid_m.centroid_m.x_m, -0.5);
  EXPECT_DOUBLE_EQ(nonzero_moments.liquid_m.centroid_m.y_m, 1.0);
  EXPECT_DOUBLE_EQ(nonzero_moments.liquid_m.centroid_m.z_m, -8.0);
  EXPECT_DOUBLE_EQ(nonzero_moments.gas_m.volume_m, 0.80);
  EXPECT_DOUBLE_EQ(nonzero_moments.gas_m.centroid_m.x_m, 0.5);
  EXPECT_DOUBLE_EQ(nonzero_moments.gas_m.centroid_m.y_m, -1.0);
  EXPECT_DOUBLE_EQ(nonzero_moments.gas_m.centroid_m.z_m, 8.0);

  blank_moments += nonzero_moments;
  EXPECT_DOUBLE_EQ(blank_moments.liquid_m.volume_m, 0.25);
  EXPECT_DOUBLE_EQ(blank_moments.liquid_m.centroid_m.x_m, -0.5);
  EXPECT_DOUBLE_EQ(blank_moments.liquid_m.centroid_m.y_m, 1.0);
  EXPECT_DOUBLE_EQ(blank_moments.liquid_m.centroid_m.z_m, -8.0);
  EXPECT_DOUBLE_EQ(blank_moments.gas_m.volume_m, 0.80);
  EXPECT_DOUBLE_EQ(blank_moments.gas_m.centroid_m.x_m, 0.5);
  EXPECT_DOUBLE_EQ(blank_moments.gas_m.centroid_m.y_m, -1.0);
  EXPECT_DOUBLE_EQ(blank_moments.gas_m.centroid_m.z_m, 8.0);
}

TEST(GeometryClasses, Reconstruction) {
  Reconstruction reconstruction;
  EXPECT_EQ(reconstruction.isFlipped(), false);
  reconstruction.flip_cut_m = -1.0;
  EXPECT_EQ(reconstruction.isFlipped(), true);
}

TEST(GeometryClasses, Tet) {
  Tet blank_tet;
  for (int v = 0; v < 4; ++v) {
    EXPECT_DOUBLE_EQ(blank_tet.vertex_m[v].x_m, 0.0);
    EXPECT_DOUBLE_EQ(blank_tet.vertex_m[v].y_m, 0.0);
    EXPECT_DOUBLE_EQ(blank_tet.vertex_m[v].z_m, 0.0);
  }

  Tet set_tet(Pt(-1.0, -1.0, -1.0), Pt(-1.0, 0.0, 0.0), Pt(0.0, 0.0, -1.0),
              Pt(0.0, 0.0, 0.0));
  double vol = set_tet.volume();
  Pt centroid = set_tet.centroid();
  EXPECT_DOUBLE_EQ(vol, 1.0 / 6.0);
  EXPECT_DOUBLE_EQ(centroid.x_m, -0.5);
  EXPECT_DOUBLE_EQ(centroid.y_m, -0.25);
  EXPECT_DOUBLE_EQ(centroid.z_m, -0.5);

  GeometricMoments tet_moments = set_tet.moments();
  EXPECT_DOUBLE_EQ(tet_moments.volume_m, vol);
  EXPECT_DOUBLE_EQ(tet_moments.centroid_m.x_m, centroid.x_m);
  EXPECT_DOUBLE_EQ(tet_moments.centroid_m.y_m, centroid.y_m);
  EXPECT_DOUBLE_EQ(tet_moments.centroid_m.z_m, centroid.z_m);
}
}  // namespace
