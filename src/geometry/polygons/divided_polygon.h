// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_DIVIDED_POLYGON_H_
#define SRC_GEOMETRY_POLYGONS_DIVIDED_POLYGON_H_

#include <algorithm>
#include <array>

#include "src/geometry/general/pt.h"
#include "src/geometry/half_edge_structures/half_edge_polygon.h"
#include "src/geometry/polygons/general_polygon.h"
#include "src/geometry/polygons/polygon_storage_types.h"
#include "src/geometry/polygons/tri.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"
#include "src/helpers/serializer.h"
#include "src/moments/volume_moments.h"

namespace IRL {

template <class MaskedType>
class MaskStripper;

/// \brief A polygon that knows its centroid and
/// can return the triangles that it is comprised of.
template <class Derived, class VertexType>
class DividedPolygonSpecialization
    : public GeneralPolygon<Derived, VertexType,
                            ProxyTri<MaskStripper<Derived>>> {
  Derived& getDerived(void) { return static_cast<Derived&>(*this); }
  const Derived& getDerived(void) const {
    return static_cast<const Derived&>(*this);
  }

 public:
  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const;

  std::array<UnsignedIndex_t, 3> getSimplexIndicesFromDecomposition(
      UnsignedIndex_t a_tri_number_to_get) const;

  /// \brief Returns the triangle.
  ProxyTri<MaskStripper<Derived>> getSimplexFromDecomposition(
      UnsignedIndex_t a_tri_number_to_get) const;

  HalfEdgePolygon<VertexType> generateHalfEdgeVersion(void) const;

  template <class HalfEdgePolygonClass>
  inline void setHalfEdgeVersion(
      HalfEdgePolygonClass* a_half_edge_version) const;

  void resetCentroid(void);

  Pt calculateCentroid(void) const;

  VolumeMoments calculateMoments(void) const;

 protected:
  VertexType& getCentroid(void);
  const VertexType& getCentroid(void) const;

 private:
  Pt calculatePolygonCentroid(void) const;
};

template <class VertexType>
class ExpandableDividedPolygon
    : public DelayedExpandableStorageAndStoredPlane<
          ExpandableDividedPolygon<VertexType>, VertexType, 1>,
      public DividedPolygonSpecialization<ExpandableDividedPolygon<VertexType>,
                                          VertexType> {
  friend DelayedExpandableStorageAndStoredPlane<
      ExpandableDividedPolygon<VertexType>, VertexType, 1>;

 public:
  using DelayedExpandableStorageAndStoredPlane<
      ExpandableDividedPolygon<VertexType>, VertexType,
      1>::DelayedExpandableStorageAndStoredPlane;

  ExpandableDividedPolygon(void) = default;

  template <class AnyDerivedClass, class AnySimplexType>
  static ExpandableDividedPolygon fromPolygon(
      const GeneralPolygon<AnyDerivedClass, VertexType, AnySimplexType>&
          a_polygon);

 private:
  template <class AnyDerivedClass, class AnySimplexType>
  ExpandableDividedPolygon(const GeneralPolygon<AnyDerivedClass, VertexType,
                                                AnySimplexType>& a_polygon);
};

///////////////////////////////////////////
// Predefined types
using DividedPolygon = ExpandableDividedPolygon<Pt>;
}  // namespace IRL

#include "src/geometry/polygons/divided_polygon.tpp"

#endif  // SRC_GEOMETRY_POLYGONS_DIVIDED_POLYGON_H_
