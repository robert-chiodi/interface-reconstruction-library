// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PT_LIST_TPP_
#define SRC_GEOMETRY_GENERAL_PT_LIST_TPP_

namespace IRL {

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
PtList<VertexType, kMaxNumberOfPts>::PtList(
    std::initializer_list<VertexType> a_list)
    : pt_list_m() {
  assert(a_list.size() == pt_list_m.size());
  std::copy(a_list.begin(), a_list.end(), pt_list_m.begin());
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
PtList<VertexType, kMaxNumberOfPts>
PtList<VertexType, kMaxNumberOfPts>::fromRawPtPointer(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts) {
  assert(a_array_of_pts != nullptr);
  return PtList(a_number_of_pts, a_array_of_pts);
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
PtList<VertexType, kMaxNumberOfPts>
PtList<VertexType, kMaxNumberOfPts>::fromRawDoublePointer(
    const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs) {
  return PtList(a_number_of_pts, a_array_of_locs);
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
constexpr UnsignedIndex_t PtList<VertexType, kMaxNumberOfPts>::getNumberOfPts(
    void) {
  return kMaxNumberOfPts;
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
const VertexType* PtList<VertexType, kMaxNumberOfPts>::getPtList(void) const {
  return pt_list_m.data();
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
VertexType& PtList<VertexType, kMaxNumberOfPts>::operator[](
    const UnsignedIndex_t a_index) {
  assert(a_index < kMaxNumberOfPts);
  return pt_list_m[a_index];
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
const VertexType& PtList<VertexType, kMaxNumberOfPts>::operator[](
    const UnsignedIndex_t a_index) const {
  assert(a_index < kMaxNumberOfPts);
  return pt_list_m[a_index];
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
IRL::Pt PtList<VertexType, kMaxNumberOfPts>::getLowerLimits(void) const {
  Pt pt_to_return(DBL_MAX, DBL_MAX, DBL_MAX);
  for (const auto& pt : pt_list_m) {
    pt_to_return[0] = std::min(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::min(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::min(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
IRL::Pt PtList<VertexType, kMaxNumberOfPts>::getUpperLimits(void) const {
  Pt pt_to_return(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (const auto& pt : pt_list_m) {
    pt_to_return[0] = std::max(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::max(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::max(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
typename PtList<VertexType, kMaxNumberOfPts>::iterator
PtList<VertexType, kMaxNumberOfPts>::begin(void) noexcept {
  return pt_list_m.begin();
}
template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
typename PtList<VertexType, kMaxNumberOfPts>::const_iterator
PtList<VertexType, kMaxNumberOfPts>::begin(void) const noexcept {
  return this->cbegin();
}
template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
typename PtList<VertexType, kMaxNumberOfPts>::const_iterator
PtList<VertexType, kMaxNumberOfPts>::cbegin(void) const noexcept {
  return pt_list_m.cbegin();
}
template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
typename PtList<VertexType, kMaxNumberOfPts>::iterator
PtList<VertexType, kMaxNumberOfPts>::end(void) noexcept {
  return pt_list_m.end();
}
template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
typename PtList<VertexType, kMaxNumberOfPts>::const_iterator
PtList<VertexType, kMaxNumberOfPts>::end(void) const noexcept {
  return this->cend();
}
template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
typename PtList<VertexType, kMaxNumberOfPts>::const_iterator
PtList<VertexType, kMaxNumberOfPts>::cend(void) const noexcept {
  return pt_list_m.cend();
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
LargeOffsetIndex_t PtList<VertexType, kMaxNumberOfPts>::getSerializedSize(
    void) const {
  return sizeof(UnsignedIndex_t) + this->getNumberOfPts() * sizeof(VertexType);
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
void PtList<VertexType, kMaxNumberOfPts>::serialize(
    ByteBuffer* a_buffer) const {
  UnsignedIndex_t number_of_points = this->getNumberOfPts();
  a_buffer->pack(&number_of_points, 1);
  for (const auto& point : pt_list_m) {
    point.serialize(a_buffer);
  }
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
void PtList<VertexType, kMaxNumberOfPts>::unpackSerialized(
    ByteBuffer* a_buffer) {
  UnsignedIndex_t number_of_points;
  a_buffer->unpack(&number_of_points, 1);
  assert(number_of_points == this->getNumberOfPts());
  for (auto& point : pt_list_m) {
    point.unpackSerialized(a_buffer);
  }
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
constexpr PtList<VertexType, kMaxNumberOfPts>::PtList(const VertexType& a_pt0,
                                                      const VertexType& a_pt1,
                                                      const VertexType& a_pt2,
                                                      const VertexType& a_pt3)
    : pt_list_m{a_pt0, a_pt1, a_pt2, a_pt3} {}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
PtList<VertexType, kMaxNumberOfPts>::PtList(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts)
    : pt_list_m() {
  assert(a_array_of_pts != nullptr);
  assert(a_number_of_pts <= kMaxNumberOfPts);
  for (UnsignedIndex_t n = 0; n < a_number_of_pts; ++n) {
    (*this)[n] = a_array_of_pts[n];
  }
}

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
PtList<VertexType, kMaxNumberOfPts>::PtList(
    const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs)
    : pt_list_m() {
  assert(a_number_of_pts <= kMaxNumberOfPts);
  for (UnsignedIndex_t n = 0; n < a_number_of_pts; ++n) {
    (*this)[n] = Pt(a_array_of_locs[3 * n], a_array_of_locs[3 * n + 1],
                    a_array_of_locs[3 * n + 2]);
  }
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_PT_LIST_TPP_
