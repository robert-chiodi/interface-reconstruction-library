// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_STORED_VERTEX_ACCESS_TPP_
#define SRC_GEOMETRY_GENERAL_STORED_VERTEX_ACCESS_TPP_

namespace IRL {

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::StoredVertexAccess(
    std::initializer_list<VertexType> a_list)
    : pt_list_m{a_list} {}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
Derived
StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::fromRawPtPointer(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts) {
  assert(a_array_of_pts != nullptr);
  return Derived(a_number_of_pts, a_array_of_pts);
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
Derived StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::
    fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                         const double* a_array_of_locs) {
  return Derived(a_number_of_pts, a_array_of_locs);
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
VertexType& StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::access(
    const UnsignedIndex_t a_index) {
  assert(a_index < kNumberOfVertices);
  return pt_list_m[a_index];
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
const VertexType&
StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::access(
    const UnsignedIndex_t a_index) const {
  assert(a_index < kNumberOfVertices);
  return pt_list_m[a_index];
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
constexpr UnsignedIndex_t StoredVertexAccess<
    Derived, VertexType, kNumberOfVertices>::getNumberOfVerticesInObject(void) {
  return kNumberOfVertices;
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
LargeOffsetIndex_t
StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::getSerializedSize(
    void) const {
  return pt_list_m.getSerializedSize();
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
void StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::serialize(
    ByteBuffer* a_buffer) const {
  pt_list_m.serialize(a_buffer);
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
void StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::
    unpackSerialized(ByteBuffer* a_buffer) {
  pt_list_m.unpackSerialized(a_buffer);
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::StoredVertexAccess(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts)
    : pt_list_m(PtList<VertexType, kNumberOfVertices>::fromRawPtPointer(
          a_number_of_pts, a_array_of_pts)) {
  assert(a_array_of_pts != nullptr);
}

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
StoredVertexAccess<Derived, VertexType, kNumberOfVertices>::StoredVertexAccess(
    const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs)
    : pt_list_m(PtList<VertexType, kNumberOfVertices>::fromRawDoublePointer(
          a_number_of_pts, a_array_of_locs)) {}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_STORED_VERTEX_ACCESS_TPP_
