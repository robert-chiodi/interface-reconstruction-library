// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_TET_H_
#define SRC_GEOMETRY_POLYHEDRONS_TET_H_

#include <float.h>

#include <cmath>
#include <initializer_list>

#include "src/geometry/general/normal.h"
#include "src/geometry/general/proxy_vertex_access.h"
#include "src/geometry/general/pt.h"
#include "src/geometry/general/pt_list.h"
#include "src/geometry/general/stored_vertex_access.h"
#include "src/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "src/geometry/polyhedrons/base_polyhedron.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/mymath.h"
#include "src/helpers/serializer.h"
#include "src/moments/volume_moments.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/planar_localizer.h"

namespace IRL {
/// \brief A tet from a collection of 4 points.

template <class GeometryType>
class ProxyTet;

template <class Derived, class VertexType>
class TetSpecialization
    : public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
 public:
  /// \brief Returns the tet.
  static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

  static constexpr std::array<UnsignedIndex_t, 4>
  getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

  ProxyTet<Derived> getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet) const;

  HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const;

  template <class HalfEdgePolyhedronType>
  inline void setHalfEdgeVersion(
      HalfEdgePolyhedronType* a_half_edge_version) const;

  /// \brief Const version of getPlanarRepresentation
  PlanarLocalizer getLocalizer(void) const;

  /// \brief Calculate and return signed volume of the tet.
  /// See Owkes & Desjardins, JCP, 2014.
  inline Volume calculateAbsoluteVolume(void) const;

  /// \brief Calculate and return volume of the tet.
  inline Volume calculateVolume(void) const;

  /// \brief Calculate sign for the tet
  /// See Owkes & Desjardins, JCP, 2014.
  inline double calculateSign(void) const;

  /// \brief Calculate and return centroid of the tet.
  inline Pt calculateCentroid(void) const;

  /// \brief Calculate and return volume weighted VolumeMoments.
  inline VolumeMoments calculateMoments() const;

 private:
};

template <class GeometryType>
class ProxyTet
    : public ProxyVertexAccess<ProxyTet<GeometryType>, GeometryType, 4>,
      public TetSpecialization<ProxyTet<GeometryType>,
                               typename GeometryType::pt_type> {
  friend ProxyVertexAccess<ProxyTet<GeometryType>, GeometryType, 4>;
 public:
  using ProxyVertexAccess<ProxyTet<GeometryType>, GeometryType,
                          4>::ProxyVertexAccess;
};

template <class VertexType>
class StoredTet
    : public StoredVertexAccess<StoredTet<VertexType>, VertexType, 4>,
      public TetSpecialization<StoredTet<VertexType>, VertexType> {
	friend StoredVertexAccess<StoredTet<VertexType>, VertexType, 4>;
 public:
  using StoredVertexAccess<StoredTet<VertexType>, VertexType,
                           4>::StoredVertexAccess;

  StoredTet(void) = default;

  template <class GeometryType>
  explicit StoredTet(const ProxyTet<GeometryType>& a_tet_proxy);

  template <class GeometryType>
  StoredTet& operator=(const ProxyTet<GeometryType>& a_tet_proxy);
};

// Predefined types
using Tet = StoredTet<Pt>;

}  // namespace IRL

#include "src/geometry/polyhedrons/tet.tpp"

#endif  // SRC_GEOMETRY_POLYHEDRONS_TET_H_
