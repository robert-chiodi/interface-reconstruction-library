// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYGONS_POLYGON_H_
#define IRL_GEOMETRY_POLYGONS_POLYGON_H_

#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/half_edge_structures/half_edge_polygon.h"
#include "irl/geometry/polygons/general_polygon.h"
#include "irl/geometry/polygons/polygon_storage_types.h"
#include "irl/geometry/polygons/stored_plane_of_existence.h"
#include "irl/geometry/polygons/tri.h"
#include "irl/helpers/byte_buffer.h"
#include "irl/helpers/geometric_cutting_helpers.h"
#include "irl/helpers/helper.h"
#include "irl/helpers/mymath.h"
#include "irl/helpers/serializer.h"
#include "irl/moments/volume_moments.h"
#include "irl/parameters/defined_types.h"
#include "irl/planar_reconstruction/planar_localizer.h"

namespace IRL {

/// \brief A polygon from a given plane.
template <class Derived, class VertexType>
class PolygonSpecialization
    : public GeneralPolygon<Derived, VertexType, ProxyTri<Derived>> {
 public:
  using pt_type = VertexType;

  HalfEdgePolygon<pt_type> generateHalfEdgeVersion(void) const;

  template <class HalfEdgePolygonType>
  inline void setHalfEdgeVersion(
      HalfEdgePolygonType* a_half_edge_version) const;

  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const;

  std::array<UnsignedIndex_t, 3> getSimplexIndicesFromDecomposition(
      UnsignedIndex_t a_tri_number_to_get) const;

  /// \brief Returns the triangle.
  ProxyTri<Derived> getSimplexFromDecomposition(
      const UnsignedIndex_t a_tri_number_to_get) const;
};

template <class VertexType>
class ExpandablePolygon
    : public ExpandableStorageAndStoredPlane<ExpandablePolygon<VertexType>,
                                             VertexType>,
      public PolygonSpecialization<ExpandablePolygon<VertexType>, VertexType> {
  friend ExpandableStorageAndStoredPlane<ExpandablePolygon<VertexType>,
                                         VertexType>;

 public:
  using ExpandableStorageAndStoredPlane<
      ExpandablePolygon<VertexType>,
      VertexType>::ExpandableStorageAndStoredPlane;

  ExpandablePolygon(void) = default;
};

///////////////////////////////////////////
// Predefined types
using Polygon = ExpandablePolygon<Pt>;
//////////////////////////////////////////

}  // namespace IRL

#include "irl/geometry/polygons/polygon.tpp"

#endif // IRL_GEOMETRY_POLYGONS_POLYGON_H_
