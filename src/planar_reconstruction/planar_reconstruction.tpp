// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_TPP_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_TPP_

namespace IRL {

template <UnsignedIndex_t kStackPlanes>
PlanarReconstruction<kStackPlanes>
PlanarReconstruction<kStackPlanes>::fromOnePlane(const Plane& a_plane) {
  return PlanarReconstruction(a_plane);
}

template <UnsignedIndex_t kStackPlanes>
PlanarReconstruction<kStackPlanes>
PlanarReconstruction<kStackPlanes>::fromTwoPlanes(const Plane& a_plane_0,
                                                  const Plane& a_plane_1) {
  return PlanarReconstruction(a_plane_0, a_plane_1);
}

template <UnsignedIndex_t kStackPlanes>
UnsignedIndex_t PlanarReconstruction<kStackPlanes>::getNumberOfPlanes(
    void) const {
  return static_cast<UnsignedIndex_t>(planes_m.size());
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::reserve(
    const UnsignedIndex_t a_number_of_future_planes) {
  planes_m.reserve(a_number_of_future_planes);
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::setNumberOfPlanes(
    const UnsignedIndex_t a_number_of_future_planes) {
  planes_m.resize(a_number_of_future_planes);
  this->checkIfStaticAllocationExceeded();
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::zeroNumberOfPlanes(void) {
  this->setNumberOfPlanes(0);
}

template <UnsignedIndex_t kStackPlanes>
Plane& PlanarReconstruction<kStackPlanes>::operator[](
    const UnsignedIndex_t a_p) {
  assert(a_p < this->getNumberOfPlanes());
  return planes_m[static_cast<std::size_t>(a_p)];
}

template <UnsignedIndex_t kStackPlanes>
const Plane& PlanarReconstruction<kStackPlanes>::operator[](
    const UnsignedIndex_t a_p) const {
  assert(a_p < this->getNumberOfPlanes());
  return planes_m[static_cast<std::size_t>(a_p)];
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::addPlane(const Plane& a_plane) {
  planes_m.push_back(a_plane);
  this->checkIfStaticAllocationExceeded();
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::addBeginningPlane(
    const Plane& a_plane) {
  planes_m.insert(planes_m.begin(), a_plane);
  this->checkIfStaticAllocationExceeded();
}

template <UnsignedIndex_t kStackPlanes>
template <class ArrayType>
void PlanarReconstruction<kStackPlanes>::setDistances(
    const ArrayType& a_distances) {
  assert(this->getNumberOfPlanes() == a_distances.size());
  for (UnsignedIndex_t p = 0; p < this->getNumberOfPlanes(); ++p) {
    planes_m[p].distance() = a_distances[p];
  }
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::removePlane(
    const UnsignedIndex_t a_p) {
  assert(a_p < this->getNumberOfPlanes());
  planes_m.erase(planes_m.begin() + a_p);
}

template <UnsignedIndex_t kStackPlanes>
LargeOffsetIndex_t PlanarReconstruction<kStackPlanes>::getSerializedSize(
    void) const {
  LargeOffsetIndex_t mysize =
      sizeof(UnsignedIndex_t);  // Size of number of planes
  for (const auto& plane : planes_m) {
    mysize += plane.getSerializedSize();
  }
  return mysize;
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::serialize(ByteBuffer* a_buffer) const {
  UnsignedIndex_t number_of_planes = this->getNumberOfPlanes();
  a_buffer->pack(&number_of_planes, 1);
  for (const auto& plane : planes_m) {
    plane.serialize(a_buffer);
  }
}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::unpackSerialized(
    ByteBuffer* a_buffer) {
  UnsignedIndex_t number_of_planes;
  a_buffer->unpack(&number_of_planes, 1);
  this->setNumberOfPlanes(number_of_planes);
  this->checkIfStaticAllocationExceeded();
  for (auto& plane : planes_m) {
    plane.unpackSerialized(a_buffer);
  }
}

template <UnsignedIndex_t kStackPlanes>
typename PlanarReconstruction<kStackPlanes>::iterator
PlanarReconstruction<kStackPlanes>::begin(void) noexcept {
  return planes_m.begin();
}
template <UnsignedIndex_t kStackPlanes>
typename PlanarReconstruction<kStackPlanes>::const_iterator
PlanarReconstruction<kStackPlanes>::begin(void) const noexcept {
  return this->cbegin();
}
template <UnsignedIndex_t kStackPlanes>
typename PlanarReconstruction<kStackPlanes>::const_iterator
PlanarReconstruction<kStackPlanes>::cbegin(void) const noexcept {
  return planes_m.cbegin();
}
template <UnsignedIndex_t kStackPlanes>
typename PlanarReconstruction<kStackPlanes>::iterator
PlanarReconstruction<kStackPlanes>::end(void) noexcept {
  return planes_m.end();
}
template <UnsignedIndex_t kStackPlanes>
typename PlanarReconstruction<kStackPlanes>::const_iterator
PlanarReconstruction<kStackPlanes>::end(void) const noexcept {
  return this->cend();
}
template <UnsignedIndex_t kStackPlanes>
typename PlanarReconstruction<kStackPlanes>::const_iterator
PlanarReconstruction<kStackPlanes>::cend(void) const noexcept {
  return planes_m.cend();
}

template <UnsignedIndex_t kStackPlanes>
PlanarReconstruction<kStackPlanes>::PlanarReconstruction(const Plane& a_plane)
    : planes_m({a_plane}) {}

template <UnsignedIndex_t kStackPlanes>
PlanarReconstruction<kStackPlanes>::PlanarReconstruction(const Plane& a_plane_0,
                                                         const Plane& a_plane_1)
    : planes_m({a_plane_0, a_plane_1}) {}

template <UnsignedIndex_t kStackPlanes>
void PlanarReconstruction<kStackPlanes>::checkIfStaticAllocationExceeded(
    void) const {
#ifndef NDEBUG_PERF
  if ((planes_m.capacity() > kStackPlanes)) {
    std::cout << "Static allocation size for SmallVector exceeded in "
                 "PlanarReconstruction. Expect performance "
                 "penalty if this happens frequently."
              << std::endl;
  }
#endif
}

template <UnsignedIndex_t kStackPlanes>
inline std::ostream& operator<<(
    std::ostream& out,
    const PlanarReconstruction<kStackPlanes>& a_reconstruction) {
  out << "Reconstruction consists of " << a_reconstruction.getNumberOfPlanes()
      << " planes\n";
  for (UnsignedIndex_t n = 0; n < a_reconstruction.getNumberOfPlanes(); ++n) {
    out << "Plane " << n << " : " << a_reconstruction[n];
  }
  return out;
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_TPP_
