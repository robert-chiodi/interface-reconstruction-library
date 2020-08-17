// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_DELAYED_EXPANDABLE_VERTEX_ACCESS_TPP_
#define SRC_GEOMETRY_GENERAL_DELAYED_EXPANDABLE_VERTEX_ACCESS_TPP_

namespace IRL {

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
DelayedExpandableVertexAccess<
    Derived, VertexType, kMaskedVertices>::DelayedExpandableVertexAccess(void)
    : pt_list_m() {
  // Add in the masked vertices
  this->insertMaskedVertices();
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    DelayedExpandableVertexAccess(std::initializer_list<VertexType> a_list)
    : pt_list_m{a_list} {
  // Add in the masked vertices
  this->insertMaskedVertices();
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
Derived DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                     const VertexType* a_array_of_pts) {
  assert(a_array_of_pts != nullptr);
  return Derived(a_number_of_pts, a_array_of_pts);
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
Derived DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                         const double* a_array_of_locs) {
  return Derived(a_number_of_pts, a_array_of_locs);
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
VertexType&
DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::access(
    const UnsignedIndex_t a_index) {
  return pt_list_m[a_index + kMaskedVertices];
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
const VertexType&
DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::access(
    const UnsignedIndex_t a_index) const {
  return pt_list_m[a_index + kMaskedVertices];
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
VertexType& DelayedExpandableVertexAccess<
    Derived, VertexType, kMaskedVertices>::unmaskedAccess(const UnsignedIndex_t
                                                              a_index) {
  return pt_list_m[a_index];
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
const VertexType& DelayedExpandableVertexAccess<
    Derived, VertexType, kMaskedVertices>::unmaskedAccess(const UnsignedIndex_t
                                                              a_index) const {
  return pt_list_m[a_index];
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
UnsignedIndex_t DelayedExpandableVertexAccess<
    Derived, VertexType, kMaskedVertices>::getNumberOfVerticesInObject(void)
    const {
  return pt_list_m.getNumberOfPts() - kMaskedVertices;
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    setNumberOfVerticesInObject(const UnsignedIndex_t a_size) {
  pt_list_m.setNumberOfPts(a_size + kMaskedVertices);
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    addVertex(const VertexType& a_vertex) {
  pt_list_m.addPt(a_vertex);
}
template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    addVertexAtIndex(const UnsignedIndex_t a_index, const VertexType& a_pt) {
  pt_list_m.addPtAtIndex(a_index + kMaskedVertices, a_pt);
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    reserve(const UnsignedIndex_t a_size) {
  pt_list_m.reserve(a_size + kMaskedVertices);
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType,
                                   kMaskedVertices>::removeLastVertex(void) {
  pt_list_m.removeLastPt();
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    removePt(const UnsignedIndex_t a_index) {
  assert(a_index < this->getNumberOfVerticesInObject());
  pt_list_m.removePt(a_index + kMaskedVertices);
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
LargeOffsetIndex_t DelayedExpandableVertexAccess<
    Derived, VertexType, kMaskedVertices>::getSerializedSize(void) const {
  return pt_list_m.getSerializedSize();
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    serialize(ByteBuffer* a_buffer) const {
  pt_list_m.serialize(a_buffer);
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    unpackSerialized(ByteBuffer* a_buffer) {
  pt_list_m.unpackSerialized(a_buffer);
  // Add in the masked vertices
  this->insertMaskedVertices();
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    DelayedExpandableVertexAccess(const UnsignedIndex_t a_number_of_pts,
                                  const VertexType* a_array_of_pts)
    : pt_list_m(ExpandablePtList<VertexType>::fromRawPtPointer(
          a_number_of_pts, a_array_of_pts)) {
  assert(a_array_of_pts != nullptr);
  // Add in the masked vertices
  this->insertMaskedVertices();
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
DelayedExpandableVertexAccess<Derived, VertexType, kMaskedVertices>::
    DelayedExpandableVertexAccess(const UnsignedIndex_t a_number_of_pts,
                                  const double* a_array_of_locs)
    : pt_list_m(ExpandablePtList<VertexType>::fromRawDoublePointer(
          a_number_of_pts, a_array_of_locs)) {
  // Add in the masked vertices
  this->insertMaskedVertices();
}

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
void DelayedExpandableVertexAccess<
    Derived, VertexType, kMaskedVertices>::insertMaskedVertices(void) {
  for (UnsignedIndex_t n = 0; n < kMaskedVertices; ++n) {
    pt_list_m.addPtAtIndex(n, VertexType());
  }
}

template <class MaskedType>
typename MaskStripper<MaskedType>::pt_type& MaskStripper<MaskedType>::
operator[](const UnsignedIndex_t a_index) {
  return MaskedType::unmaskedAccess(a_index);
}

template <class MaskedType>
const typename MaskStripper<MaskedType>::pt_type& MaskStripper<MaskedType>::
operator[](const UnsignedIndex_t a_index) const {
  return MaskedType::unmaskedAccess(a_index);
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_DELAYED_EXPANDABLE_VERTEX_ACCESS_TPP_
