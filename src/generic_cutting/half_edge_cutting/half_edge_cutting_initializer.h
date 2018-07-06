// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_H_
#define SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_H_

/// \brief Decomposing of a half-edge data structure representing an
/// initial object into one separated by a series of planar reconstructions.
///
/// This cutting takes an original object, such as Tetrahedron or Dodecahedron
/// and applies its setHalfEdgeVersion method to on a static allocation of a
/// HalfEdgePolygon or HalfEdgePolyhedron. This is done through lazy evaluation
/// in the function setHalfEdgeStructure, since the memory addresses for each
/// HalfEdge/VertexLocation/Vertex/Face will be the same, with only the
/// VertexLocation's Pt changing. A SegmentedHalfEdgePolygon or
/// SegmentedHalfEdgePolyhedron is then generated from this, which is just
/// a collection of vertices and faces from this over-arching "complete"
/// HalfEdgePolygon/HalfEdgePolyhedron. This is then cut by the
/// planar reconstruction objects through the call to
/// getVolumeMoments(...), where a pointer to the Segmented and
/// complete polytope are provided, along with the reconstruction.
/// If only the area underneath the reconstruction is needed, the
/// SegmentedPolytope is changed in place, where the Faces and vertices are
/// changed, with new ones added, representing the new half-edge structure of
/// the object that lays completely underneath the reconstructions. If instead
/// the objects on both sides of the plane are needed, such as during linked
/// cutting where volumes are distributed across a system of reconstructions
/// with known representations, then new SegmentedHalfEdgePolytope objects are
/// generated to represent the clipped (above plane) portions. During the
/// cutting of these half edge structures, all new
/// HalfEdge/VertexLocation/Vertex/Face objects are taken from the
/// allocation of the original HalfEdgePolytope, and therefore this
/// must have stable memory so that pointers remain valid during cutting.
/// It is for this reason that the original HalfEdgePolytope is referred
/// to as "complete", since at the end it will contain the complete polytope,
/// subdivided amonst the reconstructions, an accumulation of the
/// Segmented versions generated during cutting. Note, since Segmented versions
/// essentially store all data on the original HalfEdgePolytope, with just
/// pointers to the data it owns there, the Segmented versions can be made
/// invalid if the complete_polytope is modified directly, so it's best
/// not to do this if SegmentedPolytopes you hope to use still exist.

#include "src/generic_cutting/general/class_classifications.h"
#include "src/geometry/general/geometry_type_traits.h"
#include "src/helpers/SFINAE_boiler_plate.h"

namespace IRL {

// Overarching function called from generic_cutting.h
template <class ReturnType, class EncompassingType, class ReconstructionType>
inline ReturnType cutThroughHalfEdgeStructures(
    const EncompassingType& a_polytope,
    const ReconstructionType& a_reconstruction);

}  // namespace IRL

#include "src/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"

#endif  // SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_H_
