// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_POLYGON_STORAGE_TYPES_TPP_
#define SRC_GEOMETRY_POLYGONS_POLYGON_STORAGE_TYPES_TPP_

namespace IRL {

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
PolygonConstructors<Derived, VertexType, StorageType, PlaneOfExistenceType>::
    PolygonConstructors(std::initializer_list<VertexType> a_list,
                        const Plane& a_plane_of_existence)
    : StorageType{a_list}, PlaneOfExistenceType(a_plane_of_existence) {}

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
Derived
PolygonConstructors<Derived, VertexType, StorageType, PlaneOfExistenceType>::
    fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                     const VertexType* a_array_of_pts) {
  assert(a_array_of_pts != nullptr);
  return Derived(a_number_of_pts, a_array_of_pts);
}

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
Derived
PolygonConstructors<Derived, VertexType, StorageType, PlaneOfExistenceType>::
    fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                         const double* a_array_of_locs) {
  return Derived(a_number_of_pts, a_array_of_locs);
}

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
Derived
PolygonConstructors<Derived, VertexType, StorageType, PlaneOfExistenceType>::
    fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                     const VertexType* a_array_of_pts,
                     const Plane& a_plane_of_existence) {
  assert(a_array_of_pts != nullptr);
  return Derived(a_number_of_pts, a_array_of_pts, a_plane_of_existence);
}

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
Derived
PolygonConstructors<Derived, VertexType, StorageType, PlaneOfExistenceType>::
    fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                         const double* a_array_of_locs,
                         const Plane& a_plane_of_existence) {
  return Derived(a_number_of_pts, a_array_of_locs, a_plane_of_existence);
}

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
PolygonConstructors<Derived, VertexType, StorageType, PlaneOfExistenceType>::
    PolygonConstructors(const UnsignedIndex_t a_number_of_pts,
                        const VertexType* a_array_of_pts,
                        const Plane& a_plane_of_existence)
    : StorageType(a_number_of_pts, a_array_of_pts),
      PlaneOfExistenceType(a_plane_of_existence) {}

template <class Derived, class VertexType, class StorageType,
          class PlaneOfExistenceType>
PolygonConstructors<Derived, VertexType, StorageType, PlaneOfExistenceType>::
    PolygonConstructors(const UnsignedIndex_t a_number_of_pts,
                        const double* a_array_of_locs,
                        const Plane& a_plane_of_existence)
    : StorageType(a_number_of_pts, a_array_of_locs),
      PlaneOfExistenceType(a_plane_of_existence) {}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                         PlaneOfExistenceType>::ProxyPolygonConstructors(void)
    : ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>(),
      PlaneOfExistenceType() {}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                         PlaneOfExistenceType>::
    ProxyPolygonConstructors(const GeometryType& a_base_geometry,
                             std::initializer_list<UnsignedIndex_t> a_list)
    : ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>(
          a_base_geometry, a_list),
      PlaneOfExistenceType(a_base_geometry.getPlaneOfExistence()) {}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                         PlaneOfExistenceType>::
    ProxyPolygonConstructors(
        const GeometryType& a_base_geometry,
        std::array<UnsignedIndex_t, kNumberOfVertices> a_list)
    : ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>(
          a_base_geometry, a_list),
      PlaneOfExistenceType(a_base_geometry.getPlaneOfExistence()) {}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                         PlaneOfExistenceType>::
    ProxyPolygonConstructors(const GeometryType& a_base_geometry,
                             std::initializer_list<UnsignedIndex_t> a_list,
                             const Plane& a_plane_of_existence)
    : ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>(
          a_base_geometry, a_list),
      PlaneOfExistenceType(a_plane_of_existence) {}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                         PlaneOfExistenceType>::
    ProxyPolygonConstructors(
        const GeometryType& a_base_geometry,
        std::array<UnsignedIndex_t, kNumberOfVertices> a_list,
        const Plane& a_plane_of_existence)
    : ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>(
          a_base_geometry, a_list),
      PlaneOfExistenceType(a_plane_of_existence) {}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
Derived ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                                 PlaneOfExistenceType>::
    fromNoExistencePlane(const GeometryType& a_base_geometry,
                         std::initializer_list<UnsignedIndex_t> a_list) {
  return Derived(a_base_geometry, a_list, nullptr);
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
Derived ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                                 PlaneOfExistenceType>::
    fromNoExistencePlane(
        const GeometryType& a_base_geometry,
        std::array<UnsignedIndex_t, kNumberOfVertices> a_list) {
  return Derived(a_base_geometry, a_list, nullptr);
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                         PlaneOfExistenceType>::
    ProxyPolygonConstructors(const GeometryType& a_base_geometry,
                             std::initializer_list<UnsignedIndex_t> a_list,
                             Plane* a_plane_of_existence)
    : ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>(
          a_base_geometry, a_list),
      PlaneOfExistenceType(a_plane_of_existence) {}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices,
          class PlaneOfExistenceType>
ProxyPolygonConstructors<Derived, GeometryType, kNumberOfVertices,
                         PlaneOfExistenceType>::
    ProxyPolygonConstructors(
        const GeometryType& a_base_geometry,
        std::array<UnsignedIndex_t, kNumberOfVertices> a_list,
        Plane* a_plane_of_existence)
    : ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>(
          a_base_geometry, a_list),
      PlaneOfExistenceType(a_plane_of_existence) {}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_POLYGON_STORAGE_TYPES_TPP_
