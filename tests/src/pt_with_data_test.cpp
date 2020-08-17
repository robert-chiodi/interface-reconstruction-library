// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/pt_with_data.h"

#include <algorithm>
#include <random>

#include "gtest/gtest.h"

#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/new_pt_calculation_functors.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/moments/volume_moments_and_doubles.h"

namespace {

using namespace IRL;

TEST(PtWithData, Datastorage) {
  auto first_point =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(0.0, 0.0, 0.0), std::array<double, 3>{0.0, 0.0, 0.0});
  auto second_point =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(1.0, 1.0, 1.0), std::array<double, 3>{1.0, 1.0, 1.0});
  auto half_intersection =
      first_point.fromEdgeIntersection(first_point, -0.5, second_point, 0.5);
  EXPECT_NEAR(half_intersection[0], 0.5, 1.0e-14);
  EXPECT_NEAR(half_intersection[1], 0.5, 1.0e-14);
  EXPECT_NEAR(half_intersection[2], 0.5, 1.0e-14);
  EXPECT_NEAR(half_intersection.getData()[0], 0.5, 1.0e-14);
  EXPECT_NEAR(half_intersection.getData()[1], 0.5, 1.0e-14);
  EXPECT_NEAR(half_intersection.getData()[2], 0.5, 1.0e-14);
}

TEST(PtWithData, RecursiveSimplexCutting) {
  std::array<PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>, 4>
      tet_pt_list;

  tet_pt_list[0] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(1.0, 1.0, 0.0), std::array<double, 3>{1.0, 1.0, 1.0});
  tet_pt_list[1] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(1.0, 0.0, 1.0), std::array<double, 3>{0.0, 0.0, 0.0});
  tet_pt_list[2] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(1.0, 0.0, -1.0), std::array<double, 3>{0.0, 0.0, 0.0});
  tet_pt_list[3] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(0.0, 0.0, 0.0), std::array<double, 3>{2.0, 2.0, 2.0});

  StoredTet<PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>>
      formed_tet(
          {tet_pt_list[0], tet_pt_list[1], tet_pt_list[2], tet_pt_list[3]});
  auto cutting_plane =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, -1.0, 0.0), -0.5));

  auto volume_moments = getNormalizedVolumeMoments<VolumeMomentsAndDoubles<3>,
                                                   RecursiveSimplexCutting>(
      formed_tet, cutting_plane);

  Tet upper_tet({tet_pt_list[0].getPt(), Pt(1.0, 0.5, 0.5), Pt(1.0, 0.5, -0.5),
                 Pt(0.5, 0.5, 0.0)});
  auto upper_tet_volume = upper_tet.calculateVolume();
  auto upper_tet_centroid = upper_tet.calculateCentroid();

  EXPECT_NEAR(volume_moments.volume(), upper_tet_volume, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[0], upper_tet_centroid[0], 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[1], upper_tet_centroid[1], 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[2], upper_tet_centroid[2], 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[0], 7.0 / 8.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[1], 7.0 / 8.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[2], 7.0 / 8.0, 1.0e-14);
}

TEST(PtWithData, HalfEdgeCutting) {
  std::array<PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>, 4>
      tet_pt_list;

  tet_pt_list[0] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(1.0, 1.0, 0.0), std::array<double, 3>{1.0, 1.0, 1.0});
  tet_pt_list[1] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(1.0, 0.0, 1.0), std::array<double, 3>{0.0, 0.0, 0.0});
  tet_pt_list[2] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(1.0, 0.0, -1.0), std::array<double, 3>{0.0, 0.0, 0.0});
  tet_pt_list[3] =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>(
          Pt(0.0, 0.0, 0.0), std::array<double, 3>{2.0, 2.0, 2.0});

  StoredTet<PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>>
      formed_tet(
          {tet_pt_list[0], tet_pt_list[1], tet_pt_list[2], tet_pt_list[3]});

  auto cutting_plane =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, -1.0, 0.0), -0.5));

  auto volume_moments =
      getNormalizedVolumeMoments<VolumeMomentsAndDoubles<3>, HalfEdgeCutting>(
          formed_tet, cutting_plane);

  Tet upper_tet({tet_pt_list[0].getPt(), Pt(1.0, 0.5, 0.5), Pt(1.0, 0.5, -0.5),
                 Pt(0.5, 0.5, 0.0)});
  auto upper_tet_volume = upper_tet.calculateVolume();
  auto upper_tet_centroid = upper_tet.calculateCentroid();

  EXPECT_NEAR(volume_moments.volume(), upper_tet_volume, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[0], upper_tet_centroid[0], 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[1], upper_tet_centroid[1], 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[2], upper_tet_centroid[2], 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[0], 7.0 / 8.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[1], 7.0 / 8.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[2], 7.0 / 8.0, 1.0e-14);
}

TEST(PtWithData, AveragedVertex) {
  using MyPtType =
      PtWithDoublesStatelessFunctor<LinearInterpolation_Functor, 3>;
  using ArrType = std::array<double, 3>;
  std::array<ArrType, 8> rc_data_list;
  auto data_unit_cell =
      StoredRectangularCuboid<MyPtType>::fromOtherPolytope(unit_cell);
  rc_data_list[0] = ArrType{10.0, 0.0, -4.0};
  rc_data_list[1] = ArrType{10.0, 2.0, -4.0};
  rc_data_list[2] = ArrType{10.0, 2.0, 3.0};
  rc_data_list[3] = ArrType{10.0, 0.0, 3.0};
  rc_data_list[4] = ArrType{-10.0, 0.0, -4.0};
  rc_data_list[5] = ArrType{-10.0, 2.0, -4.0};
  rc_data_list[6] = ArrType{-10.0, 2.0, 3.0};
  rc_data_list[7] = ArrType{-10.0, 0.0, 3.0};
  for (UnsignedIndex_t n = 0; n < 8; ++n) {
    data_unit_cell[n].getData() = rc_data_list[n];
  }
  auto moments_and_data = data_unit_cell.calculateVolumeMomentsAndDoubles();
  moments_and_data.normalizeByVolume();
  EXPECT_NEAR(moments_and_data.volume(), 1.0, 1.0e-14);
  EXPECT_NEAR(moments_and_data.centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(moments_and_data.centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(moments_and_data.centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(moments_and_data.data()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(moments_and_data.data()[1], 1.0, 1.0e-14);
  EXPECT_NEAR(moments_and_data.data()[2], -0.5, 1.0e-14);

  auto cutting_plane =
      PlanarSeparator::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), 0.0));
  auto volume_moments = getNormalizedVolumeMoments<VolumeMomentsAndDoubles<3>>(
      data_unit_cell, cutting_plane);
  EXPECT_NEAR(volume_moments.volume(), 0.5, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[0], -0.25, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[0], -5.0, 1.0e-14);

  cutting_plane[0] = Plane(Normal(0.0, 1.0, 0.0), 0.0);
  volume_moments = getNormalizedVolumeMoments<VolumeMomentsAndDoubles<3>>(
      data_unit_cell, cutting_plane);
  EXPECT_NEAR(volume_moments.volume(), 0.5, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[1], -0.25, 1.0e-14);
  EXPECT_NEAR(volume_moments.centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments.data()[1], 0.5, 1.0e-14);

  cutting_plane =
      PlanarSeparator::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), -0.25));
  auto sep_volume_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMomentsAndDoubles<3>>>(
          data_unit_cell, cutting_plane);
  EXPECT_NEAR(sep_volume_moments[0].volume(), 0.25, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[0].centroid()[0], -0.375, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[0].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[0].centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[0].data()[0], -7.5, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[1].volume(), 0.75, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[1].centroid()[0], 0.125, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[1].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[1].centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(sep_volume_moments[1].data()[0], 2.5, 1.0e-14);
}

}  // namespace
