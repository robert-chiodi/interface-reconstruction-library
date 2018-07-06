// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CONCAVE_BOX_H_
#define SRC_GEOMETRY_POLYHEDRONS_CONCAVE_BOX_H_

#include <initializer_list>

#include "src/geometry/general/stored_vertex_access.h"
#include "src/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "src/geometry/polyhedrons/base_polyhedron.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class Derived, class VertexType>
class ConcaveBoxSpecialization
    : public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
 public:
  HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const;

  template <class HalfEdgePolyhedronType>
  void setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const;

  static constexpr std::array<UnsignedIndex_t, 4>
  getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

  static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

  ProxyTet<Derived> getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet) const;
};

template <class VertexType>
class StoredConcaveBox
    : public StoredVertexAccess<StoredConcaveBox<VertexType>, VertexType, 14>,
      public ConcaveBoxSpecialization<StoredConcaveBox<VertexType>,
                                      VertexType> {
  friend StoredVertexAccess<StoredConcaveBox<VertexType>, VertexType, 14>;

 public:
  using StoredVertexAccess<StoredConcaveBox<VertexType>, VertexType,
                           14>::StoredVertexAccess;

  StoredConcaveBox(void) = default;
};

// Predefined types
using ConcaveBox = StoredConcaveBox<Pt>;

}  // namespace IRL

#include "src/geometry/polyhedrons/concave_box.tpp"
#endif  // SRC_GEOMETRY_POLYHEDRONS_CONCAVE_BOX_H_
