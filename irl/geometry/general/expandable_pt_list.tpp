// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_EXPANDABLE_PT_LIST_TPP_
#define IRL_GEOMETRY_GENERAL_EXPANDABLE_PT_LIST_TPP_

namespace IRL {

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
ExpandablePtList<VertexType, kStaticAllocSize>
ExpandablePtList<VertexType, kStaticAllocSize>::fromRawPtPointer(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts) {
  assert(a_array_of_pts != nullptr);
  return ExpandablePtList(a_number_of_pts, a_array_of_pts);
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
ExpandablePtList<VertexType, kStaticAllocSize>
ExpandablePtList<VertexType, kStaticAllocSize>::fromRawDoublePointer(
    const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs) {
  return ExpandablePtList(a_number_of_pts, a_array_of_locs);
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
UnsignedIndex_t ExpandablePtList<VertexType, kStaticAllocSize>::getNumberOfPts(
    void) const {
  return static_cast<UnsignedIndex_t>(pt_list_m.size());
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::setNumberOfPts(
    const UnsignedIndex_t a_number) {
  pt_list_m.resize(a_number);
  this->checkIfStaticAllocationExceeded();
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::reserve(
    const UnsignedIndex_t a_number) {
  pt_list_m.reserve(a_number);
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::addPt(
    const VertexType& a_pt) {
  pt_list_m.push_back(a_pt);
  this->checkIfStaticAllocationExceeded();
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::addPtAtIndex(
    const UnsignedIndex_t a_index, const VertexType& a_pt) {
  pt_list_m.insert(this->begin() + a_index, a_pt);
  this->checkIfStaticAllocationExceeded();
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
VertexType& ExpandablePtList<VertexType, kStaticAllocSize>::operator[](
    const UnsignedIndex_t a_index) {
  assert(a_index < this->getNumberOfPts());
  return pt_list_m[a_index];
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
const VertexType& ExpandablePtList<VertexType, kStaticAllocSize>::operator[](
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->getNumberOfPts());
  return pt_list_m[a_index];
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::removePt(
    const UnsignedIndex_t a_index) {
  assert(a_index < this->getNumberOfPts());
  pt_list_m.erase(pt_list_m.begin() + a_index);
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::removeLastPt(void) {
  pt_list_m.pop_back();
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
VertexType ExpandablePtList<VertexType, kStaticAllocSize>::getLowerLimits(
    void) const {
  VertexType pt_to_return(DBL_MAX, DBL_MAX, DBL_MAX);
  for (const auto& pt : pt_list_m) {
    pt_to_return[0] = std::min(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::min(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::min(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
VertexType ExpandablePtList<VertexType, kStaticAllocSize>::getUpperLimits(
    void) const {
  VertexType pt_to_return(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (const auto& pt : pt_list_m) {
    pt_to_return[0] = std::max(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::max(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::max(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
typename ExpandablePtList<VertexType, kStaticAllocSize>::iterator
ExpandablePtList<VertexType, kStaticAllocSize>::begin(void) noexcept {
  return pt_list_m.begin();
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
typename ExpandablePtList<VertexType, kStaticAllocSize>::const_iterator
ExpandablePtList<VertexType, kStaticAllocSize>::begin(void) const noexcept {
  return this->cbegin();
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
typename ExpandablePtList<VertexType, kStaticAllocSize>::const_iterator
ExpandablePtList<VertexType, kStaticAllocSize>::cbegin(void) const noexcept {
  return pt_list_m.cbegin();
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
typename ExpandablePtList<VertexType, kStaticAllocSize>::iterator
ExpandablePtList<VertexType, kStaticAllocSize>::end(void) noexcept {
  return pt_list_m.end();
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
typename ExpandablePtList<VertexType, kStaticAllocSize>::const_iterator
ExpandablePtList<VertexType, kStaticAllocSize>::end(void) const noexcept {
  return this->cend();
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
typename ExpandablePtList<VertexType, kStaticAllocSize>::const_iterator
ExpandablePtList<VertexType, kStaticAllocSize>::cend(void) const noexcept {
  return pt_list_m.cend();
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
LargeOffsetIndex_t
ExpandablePtList<VertexType, kStaticAllocSize>::getSerializedSize(void) const {
  LargeOffsetIndex_t mysize =
      sizeof(UnsignedIndex_t);  // Size of number of points
  for (const auto& point : pt_list_m) {
    mysize += point.getSerializedSize();
  }
  return mysize;
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::serialize(
    ByteBuffer* a_buffer) const {
  UnsignedIndex_t number_of_points = this->getNumberOfPts();
  a_buffer->pack(&number_of_points, 1);
  for (const auto& point : pt_list_m) {
    point.serialize(a_buffer);
  }
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<VertexType, kStaticAllocSize>::unpackSerialized(
    ByteBuffer* a_buffer) {
  UnsignedIndex_t number_of_points;
  a_buffer->unpack(&number_of_points, 1);
  this->setNumberOfPts(number_of_points);
  this->checkIfStaticAllocationExceeded();
  for (auto& point : pt_list_m) {
    point.unpackSerialized(a_buffer);
  }
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
ExpandablePtList<VertexType, kStaticAllocSize>::ExpandablePtList(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts)
    : pt_list_m(a_number_of_pts) {
  assert(a_array_of_pts != nullptr);
  this->checkIfStaticAllocationExceeded();
  for (UnsignedIndex_t n = 0; n < a_number_of_pts; ++n) {
    (*this)[n] = a_array_of_pts[n];
  }
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
ExpandablePtList<VertexType, kStaticAllocSize>::ExpandablePtList(
    const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs)
    : pt_list_m(a_number_of_pts) {
  this->checkIfStaticAllocationExceeded();
  for (UnsignedIndex_t n = 0; n < a_number_of_pts; ++n) {
    (*this)[n] = Pt(a_array_of_locs[3 * n], a_array_of_locs[3 * n + 1],
                    a_array_of_locs[3 * n + 2]);
  }
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void ExpandablePtList<
    VertexType, kStaticAllocSize>::checkIfStaticAllocationExceeded(void) const {
#ifndef NDEBUG_PERF
  if (pt_list_m.capacity() >
      ExpandablePtList<VertexType, kStaticAllocSize>::on_stack_size) {
    std::cout << "Static allocation size for SmallVector exceeded in "
                 "ExpandablePtList. Expect performance penalty if this "
                 "happens frequently."
              << std::endl;
  }
#endif
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_EXPANDABLE_PT_LIST_TPP_
