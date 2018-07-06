// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_POLYGON_STORAGE_TYPES_H_
#define SRC_GEOMETRY_POLYGONS_POLYGON_STORAGE_TYPES_H_

#include <cassert>
#include <initializer_list>

#include "src/geometry/general/delayed_expandable_vertex_access.h"
#include "src/geometry/general/proxy_vertex_access.h"
#include "src/geometry/general/stored_vertex_access.h"
#include "src/geometry/polygons/referred_to_plane_of_existence.h"
#include "src/geometry/polygons/stored_plane_of_existence.h"
#include "src/parameters/defined_types.h"

namespace IRL {

// This is simply a way to redefine constructors needed to set both the
// vertices and the plane of existence in a way that can then just
// be inherited by polygon types that'll use these. Prevents need
// to keep rewriting the constructors each time a new polygon type is created
// that is using the same manner of vertex and plane_of_existence storage.

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
class PolygonConstructors : public StorageType, public PlaneOfExistenceType {
 public:
  using StorageType::StorageType;

  PolygonConstructors(void) : StorageType(), PlaneOfExistenceType() {}

  PolygonConstructors(std::initializer_list<VertexType> a_list,
                      const Plane& a_plane_of_existence);

  static Derived fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                                  const VertexType* a_array_of_pts);

  static Derived fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                                      const double* a_array_of_locs);

  static Derived fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                                  const VertexType* a_array_of_pts,
                                  const Plane& a_plane_of_existence);

  static Derived fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                                      const double* a_array_of_locs,
                                      const Plane& a_plane_of_existence);

 private:
  /// \brief Construct n-pts form array of pts.
  PolygonConstructors(const UnsignedIndex_t a_number_of_pts,
                      const VertexType* a_array_of_pts,
                      const Plane& a_plane_of_existence);
  /// \brief Construct n-pts form array of doubles.
  PolygonConstructors(const UnsignedIndex_t a_number_of_pts,
                      const double* a_array_of_locs,
                      const Plane& a_plane_of_existence);
};

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
class ProxyPolygonConstructors
    : public ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>,
      public PlaneOfExistenceType {
 public:
  ProxyPolygonConstructors(void);

  ProxyPolygonConstructors(const GeometryType& a_base_geometry,
                           std::initializer_list<UnsignedIndex_t> a_list);

  ProxyPolygonConstructors(const GeometryType& a_base_geometry,
                           std::initializer_list<UnsignedIndex_t> a_list,
                           const Plane& a_plane_of_existence);

  ProxyPolygonConstructors(
      const GeometryType& a_base_geometry,
      std::array<UnsignedIndex_t, kNumberOfVertices> a_list);

  ProxyPolygonConstructors(
      const GeometryType& a_base_geometry,
      std::array<UnsignedIndex_t, kNumberOfVertices> a_list,
      const Plane& a_plane_of_existence);

  static Derived fromNoExistencePlane(
      const GeometryType& a_base_geometry,
      std::initializer_list<UnsignedIndex_t> a_list);

  static Derived fromNoExistencePlane(
      const GeometryType& a_base_geometry,
      std::array<UnsignedIndex_t, kNumberOfVertices> a_list);

  // TODO(robertchiodi): This constructor should be made private, however Intel
  // compilers have trouble with finding this in ProxyTri, even when ProxyTri
  // friends this class. Need to figure out why.
  ProxyPolygonConstructors(const GeometryType& a_base_geometry,
                           std::initializer_list<UnsignedIndex_t> a_list,
                           Plane* a_plane_of_existence);

  // TODO(robertchiodi): This constructor should be made private, however Intel
  // compilers have trouble with finding this in ProxyTri, even when ProxyTri
  // friends this class. Need to figure out why.
  ProxyPolygonConstructors(
      const GeometryType& a_base_geometry,
      std::array<UnsignedIndex_t, kNumberOfVertices> a_list,
      Plane* a_plane_of_existence);

 private:
};

template <class Derived, class VertexType>
using ExpandableStorageAndStoredPlane =
    PolygonConstructors<Derived, VertexType,
                        DelayedExpandableVertexAccess<Derived, VertexType, 0>,
                        StoredPlaneOfExistence>;

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
using DelayedExpandableStorageAndStoredPlane = PolygonConstructors<
    Derived, VertexType,
    DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>,
    StoredPlaneOfExistence>;

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
using StaticStorageAndStoredPlane = PolygonConstructors<
    Derived, VertexType,
    StoredVertexAccess<Derived, VertexType, kNumberOfVertices>,
    StoredPlaneOfExistence>;

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
using ProxyStaticStorageAndRefferredToPlane =
    ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                             ReferredToPlaneOFExistence>;

}  // namespace IRL

#include "src/geometry/polygons/polygon_storage_types.tpp"

#endif  // SRC_GEOMETRY_POLYGONS_POLYGON_STORAGE_TYPES_H_
