// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_TPP_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_TPP_

namespace IRL {

inline PlanarSeparator::PlanarSeparator(void)
    : reconstruction_m(), flip_cut_m(1.0) {}

inline PlanarSeparator PlanarSeparator::fromOnePlane(const Plane& a_plane) {
  return PlanarSeparator(a_plane);
}

inline PlanarSeparator PlanarSeparator::fromTwoPlanes(
    const Plane& a_plane_0, const Plane& a_plane_1,
    const double a_flip_indicator) {
  return PlanarSeparator(a_plane_0, a_plane_1, a_flip_indicator);
}

inline UnsignedIndex_t PlanarSeparator::getNumberOfPlanes(void) const {
  return reconstruction_m.getNumberOfPlanes();
}

inline void PlanarSeparator::setNumberOfPlanes(
    const UnsignedIndex_t a_number_of_future_planes) {
  reconstruction_m.setNumberOfPlanes(a_number_of_future_planes);
}

inline void PlanarSeparator::zeroNumberOfPlanes(void) {
  reconstruction_m.zeroNumberOfPlanes();
}

inline Plane& PlanarSeparator::operator[](const UnsignedIndex_t a_p) {
  assert(a_p < this->getNumberOfPlanes());
  return reconstruction_m[a_p];
}

inline const Plane& PlanarSeparator::operator[](
    const UnsignedIndex_t a_p) const {
  assert(a_p < this->getNumberOfPlanes());
  return reconstruction_m[a_p];
}

inline void PlanarSeparator::addPlane(const Plane& a_plane) {
  reconstruction_m.addPlane(a_plane);
}

template <class ArrayType>
inline void PlanarSeparator::setDistances(
    const ArrayType& a_distances) {
  assert(this->getNumberOfPlanes() == a_distances.size());
  reconstruction_m.setDistances(a_distances);
}

inline void PlanarSeparator::removePlane(const UnsignedIndex_t a_p) {
  assert(a_p < this->getNumberOfPlanes());
  reconstruction_m.removePlane(a_p);
}

inline double PlanarSeparator::flip(void) const {
  assert(flip_cut_m == 1.0 || flip_cut_m == -1.0);
  return flip_cut_m;
}

inline void PlanarSeparator::setFlip(const double a_flip_value) {
  assert(a_flip_value == -1.0 || a_flip_value == 1.0);
  flip_cut_m = a_flip_value;
}

inline void PlanarSeparator::doNotFlipCutting(void) { flip_cut_m = 1.0; }

inline void PlanarSeparator::flipCutting(void) { flip_cut_m = -1.0; }

inline bool PlanarSeparator::isFlipped(void) const {
  return this->flip() < 0.0;
}

inline bool PlanarSeparator::isNotFlipped(void) const {
  return !this->isFlipped();
}

inline void PlanarSeparator::zeroPlanes(void) {
  this->zeroNumberOfPlanes();
  this->doNotFlipCutting();
}

inline PlanarSeparator&
PlanarSeparator::getCurrentReconstruction(void) {
  return *this;
}

inline const PlanarSeparator&
PlanarSeparator::PlanarSeparator::getCurrentReconstruction(void) const {
  return *this;
}

inline constexpr NullReconstruction PlanarSeparator::getNextReconstruction(
    void) {
  return NullReconstruction();
}

inline LargeOffsetIndex_t PlanarSeparator::getSerializedSize(void) const {
  return reconstruction_m.getSerializedSize() + sizeof(double);
}

inline void PlanarSeparator::serialize(ByteBuffer* a_buffer) const {
  reconstruction_m.serialize(a_buffer);
  a_buffer->pack(&flip_cut_m, 1);
}

inline void PlanarSeparator::unpackSerialized(ByteBuffer* a_buffer) {
  reconstruction_m.unpackSerialized(a_buffer);
  a_buffer->unpack(&flip_cut_m, 1);
}

inline PlanarSeparator::iterator PlanarSeparator::begin(void) noexcept {
  return reconstruction_m.begin();
}
inline PlanarSeparator::const_iterator PlanarSeparator::begin(void) const
    noexcept {
  return this->cbegin();
}
inline PlanarSeparator::const_iterator PlanarSeparator::cbegin(void) const
    noexcept {
  return reconstruction_m.cbegin();
}
inline PlanarSeparator::iterator PlanarSeparator::end(void) noexcept {
  return reconstruction_m.end();
}
inline PlanarSeparator::const_iterator PlanarSeparator::end(void) const
    noexcept {
  return this->cend();
}
inline PlanarSeparator::const_iterator PlanarSeparator::cend(void) const
    noexcept {
  return reconstruction_m.cend();
}

inline PlanarSeparator::PlanarSeparator(const Plane& a_plane)
    : reconstruction_m(PlanarReconstructionBase::fromOnePlane(a_plane)),
      flip_cut_m(1.0) {}

inline PlanarSeparator::PlanarSeparator(const Plane& a_plane_0,
                                        const Plane& a_plane_1,
                                        const double a_flip_indicator)
    : reconstruction_m(
          PlanarReconstructionBase::fromTwoPlanes(a_plane_0, a_plane_1)),
      flip_cut_m(a_flip_indicator) {}

inline std::ostream& operator<<(std::ostream& out,
                                const PlanarSeparator& a_reconstruction) {
  out << "PlanarSeparator: \n" << a_reconstruction.reconstruction_m;
  if (a_reconstruction.isFlipped()) {
    out << "Cutting IS flipped. \n";
  } else {
    out << "Cutting IS NOT flipped. \n";
  }
  return out;
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_TPP_
