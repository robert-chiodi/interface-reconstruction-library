// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/geometry/half_edge_structures/brep_to_half_edge.h"

#include "gtest/gtest.h"

#include "irl/geometry/polyhedrons/tet.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/helpers/geometric_cutting_helpers.h"
#include "irl/moments/volume.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"

namespace {

using namespace IRL;

TEST(BREPToHalfEdge, TetCreation) {

    std::array<std::array<UnsignedIndex_t, 3>, 4> face_mapping{{
							     {0,1,2},
							     {1,0,3},
							     {2,1,3},
							     {0,2,3}
								}};
  std::array<Pt, 4> tet_pts{{Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0),
                    Pt(0.0, 0.5, 0.0)}};
  auto made_half_edge_structure = BREPToHalfEdge<Pt>::generateHalfEdgeVersion(tet_pts, face_mapping);
  auto made_segmented = made_half_edge_structure.generateSegmentedPolyhedron();
  assert(made_segmented.checkValidHalfEdgeStructure());

  Tet tet_to_write({Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0),
                    Pt(0.0, 0.5, 0.0)});  
  auto tet_half_edge = tet_to_write.generateHalfEdgeVersion();
  auto tet_segmented_half_edge = tet_half_edge.generateSegmentedPolyhedron();
  assert(tet_segmented_half_edge.checkValidHalfEdgeStructure());

  Volume tmp_volume;
  PlanarLocalizer cutting_reconstruction =
      PlanarLocalizer::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), 0.5));
  decltype(tet_segmented_half_edge) clipped_tet;
  splitHalfEdgePolytope(&tet_segmented_half_edge, &clipped_tet, &tet_half_edge,
                        cutting_reconstruction[0]);
  decltype(tet_segmented_half_edge) clipped_made_tet;
  splitHalfEdgePolytope(&made_segmented, &clipped_made_tet, &made_half_edge_structure,
                        cutting_reconstruction[0]);
  assert(tet_segmented_half_edge.checkValidHalfEdgeStructure());
  assert(made_segmented.checkValidHalfEdgeStructure());

  EXPECT_NEAR(made_segmented.calculateVolume(), tet_segmented_half_edge.calculateVolume(),1.0e-15);
  auto made_centroid = made_segmented.calculateCentroid();
  auto tet_centroid = tet_segmented_half_edge.calculateCentroid();
  EXPECT_NEAR(made_centroid[0], tet_centroid[0],1.0e-15);
  EXPECT_NEAR(made_centroid[1], tet_centroid[1],1.0e-15);
  EXPECT_NEAR(made_centroid[2], tet_centroid[2],1.0e-15);    
}

TEST(BREPToHalfEdge, Cuboid){
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

  auto made_half_edge_structure = BREPToHalfEdge<Pt>::generateHalfEdgeVersion(cuboid_pts, face_mapping);
  auto made_segmented = made_half_edge_structure.generateSegmentedPolyhedron();
  assert(made_segmented.checkValidHalfEdgeStructure());  
  
  auto cuboid =  RectangularCuboid::fromRawPtPointer(static_cast<UnsignedIndex_t>(cuboid_pts.size()), cuboid_pts.data());
  auto cuboid_half_edge = cuboid.generateHalfEdgeVersion();
  auto cuboid_segmented_half_edge = cuboid_half_edge.generateSegmentedPolyhedron();
  assert(cuboid_segmented_half_edge.checkValidHalfEdgeStructure());  
  
  auto my_normal = Normal::normalized(1.0, 1.0, 0.0);
  PlanarLocalizer cutting_reconstruction =
    PlanarLocalizer::fromOnePlane(Plane(my_normal, my_normal*cuboid.calculateCentroid()));
  
  decltype(cuboid_segmented_half_edge) clipped_cuboid;
  splitHalfEdgePolytope(&cuboid_segmented_half_edge, &clipped_cuboid, &cuboid_half_edge,
                        cutting_reconstruction[0]);
  decltype(cuboid_segmented_half_edge) clipped_made_cuboid;
  splitHalfEdgePolytope(&made_segmented, &clipped_made_cuboid, &made_half_edge_structure,
                        cutting_reconstruction[0]);
  assert(cuboid_segmented_half_edge.checkValidHalfEdgeStructure());
  assert(made_segmented.checkValidHalfEdgeStructure());

  EXPECT_NEAR(made_segmented.calculateVolume(), cuboid_segmented_half_edge.calculateVolume(),1.0e-15);
  auto made_centroid = made_segmented.calculateCentroid();
  auto cuboid_centroid = cuboid_segmented_half_edge.calculateCentroid();
  EXPECT_NEAR(made_centroid[0], cuboid_centroid[0],1.0e-15);
  EXPECT_NEAR(made_centroid[1], cuboid_centroid[1],1.0e-15);
  EXPECT_NEAR(made_centroid[2], cuboid_centroid[2],1.0e-15);      

}
  
  
}  // namespace
