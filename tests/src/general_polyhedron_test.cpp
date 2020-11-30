// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/geometry/polyhedrons/general_polyhedron.h"

#include "irl/geometry/polyhedrons/tet.h"
#include "irl/helpers/geometric_cutting_helpers.h"
#include "irl/moments/volume.h"
#include "irl/generic_cutting/generic_cutting.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(GeneralPolyhedron, Tet) {

  std::array<std::array<UnsignedIndex_t, 3>, 4> face_mapping{{
							      {0,1,2},
							      {1,0,3},
							      {2,1,3},
							      {0,2,3}
							      }};

  PolyhedronConnectivity connectivity(face_mapping);
  std::array<Pt, 4> tet_pts{{Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0),
                    Pt(0.0, 0.5, 0.0)}};
  GeneralPolyhedron gen_poly(tet_pts, &connectivity);

  Tet tet_compile_time({Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0),
			Pt(0.0, 0.5, 0.0)});  
  PlanarLocalizer cutting_reconstruction =
      PlanarLocalizer::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), 0.5));


  auto corr_moments = getVolumeMoments<VolumeMoments, RecursiveSimplexCutting>(tet_compile_time,
								       cutting_reconstruction);
  auto gen_recursive_simplex = getVolumeMoments<VolumeMoments, RecursiveSimplexCutting>(gen_poly,
										cutting_reconstruction);
  auto gen_half_edge = getVolumeMoments<VolumeMoments, HalfEdgeCutting>(gen_poly,
									cutting_reconstruction);
 
  EXPECT_NEAR(corr_moments.volume(), gen_recursive_simplex.volume(), 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[0], gen_recursive_simplex.centroid()[0], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[1], gen_recursive_simplex.centroid()[1], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[2], gen_recursive_simplex.centroid()[2], 1.0e-15);    

  EXPECT_NEAR(corr_moments.volume(), gen_half_edge.volume(), 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[0], gen_half_edge.centroid()[0], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[1], gen_half_edge.centroid()[1], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[2], gen_half_edge.centroid()[2], 1.0e-15);

  Tet tet_compile_time_negative({Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5),
  				 Pt(0.0, 0.5, 0.0)});
  std::array<Pt, 4> negative_tet_pts{{Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5),
                    Pt(0.0, 0.5, 0.0)}};  
  auto new_gen_poly = GeneralPolyhedron(negative_tet_pts, &connectivity);

  corr_moments = getVolumeMoments<VolumeMoments, RecursiveSimplexCutting>(tet_compile_time_negative,
  								       cutting_reconstruction);
  gen_recursive_simplex = getVolumeMoments<VolumeMoments, RecursiveSimplexCutting>(new_gen_poly,
  										   cutting_reconstruction);
  gen_half_edge = getVolumeMoments<VolumeMoments, HalfEdgeCutting>(new_gen_poly,
  								   cutting_reconstruction);

  EXPECT_NEAR(corr_moments.volume(), gen_recursive_simplex.volume(), 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[0], gen_recursive_simplex.centroid()[0], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[1], gen_recursive_simplex.centroid()[1], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[2], gen_recursive_simplex.centroid()[2], 1.0e-15);
  
  EXPECT_NEAR(corr_moments.volume(), gen_half_edge.volume(), 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[0], gen_half_edge.centroid()[0], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[1], gen_half_edge.centroid()[1], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[2], gen_half_edge.centroid()[2], 1.0e-15);
  
}

TEST(GeneralPolyhedron, Cuboid){
  std::array<std::array<UnsignedIndex_t, 4>, 6> face_mapping{{
							    {0,1,2,3},
							    {1,0,4,5},
							    {2,1,5,6},
							    {3,2,6,7},
							    {0,3,7,4},
							    {7,6,5,4}
							      }};
  std::array<Pt, 8> cuboid_pts{{Pt(2.0,2.0,2.0),
				Pt(2.0,3.0,2.0),
				Pt(2.0,3.0,4.0),
				Pt(2.0,2.0,4.0),
				Pt(1.0,2.0,2.0),
				Pt(1.0,3.0,2.0),
				Pt(1.0,3.0,4.0),
				Pt(1.0,2.0,4.0)
				}};
  PolyhedronConnectivity connectivity(face_mapping); 
  GeneralPolyhedron gen_poly(cuboid_pts, &connectivity);
  
  auto cuboid =  RectangularCuboid::fromRawPtPointer(static_cast<UnsignedIndex_t>(cuboid_pts.size()), cuboid_pts.data());
  
  
  auto my_normal = Normal::normalized(1.0, 1.0, 0.0);
  PlanarLocalizer cutting_reconstruction =
    PlanarLocalizer::fromOnePlane(Plane(my_normal, my_normal*cuboid.calculateCentroid()));
  
  auto corr_moments = getVolumeMoments<VolumeMoments, RecursiveSimplexCutting>(cuboid,
									       cutting_reconstruction);
  auto gen_recursive_simplex = getVolumeMoments<VolumeMoments, RecursiveSimplexCutting>(gen_poly,
											cutting_reconstruction);
  auto gen_half_edge = getVolumeMoments<VolumeMoments, HalfEdgeCutting>(gen_poly,
									cutting_reconstruction);
  
  EXPECT_NEAR(corr_moments.volume(), gen_recursive_simplex.volume(), 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[0], gen_recursive_simplex.centroid()[0], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[1], gen_recursive_simplex.centroid()[1], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[2], gen_recursive_simplex.centroid()[2], 1.0e-15);
  
  EXPECT_NEAR(corr_moments.volume(), gen_half_edge.volume(), 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[0], gen_half_edge.centroid()[0], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[1], gen_half_edge.centroid()[1], 1.0e-15);
  EXPECT_NEAR(corr_moments.centroid()[2], gen_half_edge.centroid()[2], 1.0e-15);

}

}  // namespace
