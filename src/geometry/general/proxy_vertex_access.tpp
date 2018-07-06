// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PROXY_VERTEX_ACCESS_TPP_
#define SRC_GEOMETRY_GENERAL_PROXY_VERTEX_ACCESS_TPP_

#include <utility>

namespace IRL {

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::ProxyVertexAccess(
    const GeometryType& a_base_geometry,
    std::initializer_list<UnsignedIndex_t> a_list)
    : base_geometry_m(a_base_geometry) {
  assert(a_list.size() == index_mapping_m.size());
  std::copy(a_list.begin(), a_list.end(), index_mapping_m.begin());
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::ProxyVertexAccess(
    const GeometryType& a_base_geometry,
    std::array<UnsignedIndex_t, kNumberOfVertices> a_list)
    : base_geometry_m(a_base_geometry) {
  assert(a_list.size() == index_mapping_m.size());
  std::copy(a_list.begin(), a_list.end(), index_mapping_m.begin());
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
const typename ProxyVertexAccess<Derived, GeometryType,
                                 kNumberOfVertices>::VertexType&
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::access(
    const UnsignedIndex_t a_index) const {
  assert(a_index < index_mapping_m.size());
  return base_geometry_m[index_mapping_m[a_index]];
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
const std::array<UnsignedIndex_t, kNumberOfVertices>&
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::getIndexMapping(
    void) const {
  return index_mapping_m;
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
constexpr UnsignedIndex_t
ProxyVertexAccess<Derived, GeometryType,
                  kNumberOfVertices>::getNumberOfVerticesInObject(void) {
  return kNumberOfVertices;
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
LargeOffsetIndex_t
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::getSerializedSize(
    void) const {
  return sizeof(UnsignedIndex_t) +
         this->getNumberOfVerticesInObject() * sizeof(VertexType);
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
void ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::serialize(
    ByteBuffer* a_buffer) const {
  UnsignedIndex_t number_of_points = this->getNumberOfVerticesInObject();
  a_buffer->pack(&number_of_points, 1);
  for (UnsignedIndex_t n = 0; n < this->getNumberOfVerticesInObject(); ++n) {
    const auto& point = this->access(n);
    point.serialize(a_buffer);
  }
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::ProxyVertexAccess(
    const ProxyVertexAccess& a_other)
    : base_geometry_m(a_other.base_geometry_m) {
  index_mapping_m = a_other.index_mapping_m;
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>&
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::operator=(
    const ProxyVertexAccess& a_other) {
  assert(&base_geometry_m == &a_other.base_geometry_m);
  if (this != &a_other) {
    index_mapping_m = a_other.index_mapping_m;
  }
  return (*this);
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::ProxyVertexAccess(
    ProxyVertexAccess&& a_other)
    : base_geometry_m(a_other.base_geometry_m) {
  index_mapping_m = std::move(a_other.index_mapping_m);
}

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>&
ProxyVertexAccess<Derived, GeometryType, kNumberOfVertices>::operator=(
    ProxyVertexAccess&& a_other) {
  assert(&base_geometry_m == &a_other.base_geometry_m);
  if (this != &a_other) {
    index_mapping_m = std::move(a_other.index_mapping_m);
  }
  return (*this);
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_PROXY_VERTEX_ACCESS_TPP_
