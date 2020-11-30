// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_HEXAHEDRON_H_
#define IRL_GEOMETRY_POLYHEDRONS_HEXAHEDRON_H_

#include <float.h>

#include <cmath>
#include <initializer_list>

#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/stored_vertex_access.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "irl/geometry/polygons/tri.h"
#include "irl/geometry/polyhedrons/base_polyhedron.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/helpers/byte_buffer.h"
#include "irl/helpers/mymath.h"
#include "irl/helpers/serializer.h"
#include "irl/moments/volume.h"
#include "irl/moments/volume_moments.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/planar_localizer.h"

namespace IRL {

template <class Derived, class VertexType>
class RectangularCuboidCommon;

template <class VertexType>
class RectangularCuboidBase;

/// \brief A hexahedron class.
template <class Derived, class VertexType>
class HexahedronSpecialization
    : public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
 public:
  HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const;

  template <class HalfEdgePolyhedronType>
  inline void setHalfEdgeVersion(
      HalfEdgePolyhedronType* a_half_edge_version) const;

  /// \brief Returns a planar reconstruction that is equivalent to the
  /// HexahedronCommon.
  PlanarLocalizer getLocalizer(void) const;

  static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

  static constexpr std::array<UnsignedIndex_t, 4>
  getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

  ProxyTet<Derived> getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet) const;

 private:
};

template <class VertexType>
class StoredHexahedron
    : public StoredVertexAccess<StoredHexahedron<VertexType>, VertexType, 8>,
      public HexahedronSpecialization<StoredHexahedron<VertexType>,
                                      VertexType> {
  friend StoredVertexAccess<StoredHexahedron<VertexType>, VertexType, 8>;
 public:
  using StoredVertexAccess<StoredHexahedron<VertexType>, VertexType,
                           8>::StoredVertexAccess;

  StoredHexahedron(void) = default;
};

// Predefined types
using Hexahedron = StoredHexahedron<Pt>;
}  // namespace IRL

#include "irl/geometry/polyhedrons/hexahedron.tpp"

#endif // IRL_GEOMETRY_HEXAHEDRON_H_
