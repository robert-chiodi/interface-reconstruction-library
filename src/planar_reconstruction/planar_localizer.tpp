// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_LOCALIZER_TPP_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_LOCALIZER_TPP_

namespace IRL {

inline PlanarLocalizer PlanarLocalizer::fromOnePlane(const Plane& a_plane) {
  return PlanarLocalizer(a_plane);
}

inline PlanarLocalizer PlanarLocalizer::fromTwoPlanes(const Plane& a_plane_0,
                                                      const Plane& a_plane_1) {
  return PlanarLocalizer(a_plane_0, a_plane_1);
}

inline UnsignedIndex_t PlanarLocalizer::getNumberOfPlanes(void) const {
  return reconstruction_m.getNumberOfPlanes();
}

inline void PlanarLocalizer::setNumberOfPlanes(
    const UnsignedIndex_t a_number_of_future_planes) {
  reconstruction_m.setNumberOfPlanes(a_number_of_future_planes);
}

inline void PlanarLocalizer::zeroNumberOfPlanes(void) {
  reconstruction_m.zeroNumberOfPlanes();
}

inline Plane& PlanarLocalizer::operator[](const UnsignedIndex_t a_p) {
  assert(a_p < this->getNumberOfPlanes());
  return reconstruction_m[a_p];
}

inline const Plane& PlanarLocalizer::operator[](
    const UnsignedIndex_t a_p) const {
  assert(a_p < this->getNumberOfPlanes());
  return reconstruction_m[a_p];
}

inline void PlanarLocalizer::addPlane(const Plane& a_plane) {
  reconstruction_m.addPlane(a_plane);
}

inline void PlanarLocalizer::addBeginningPlane(const Plane& a_plane) {
  reconstruction_m.addBeginningPlane(a_plane);
}

inline void PlanarLocalizer::removePlane(const UnsignedIndex_t a_p) {
  assert(a_p < this->getNumberOfPlanes());
  reconstruction_m.removePlane(a_p);
}

inline constexpr double PlanarLocalizer::flip(void) { return 1.0; }

inline constexpr bool PlanarLocalizer::isFlipped(void) { return false; }

inline constexpr bool PlanarLocalizer::isNotFlipped(void) {
  return !PlanarLocalizer::isFlipped();
}

inline const PlanarLocalizer&
PlanarLocalizer::PlanarLocalizer::getCurrentReconstruction(void) const {
  return *this;
}

inline constexpr NullReconstruction PlanarLocalizer::getNextReconstruction(
    void) {
  return NullReconstruction();
}

inline PlanarLocalizer::iterator PlanarLocalizer::begin(void) noexcept {
  return reconstruction_m.begin();
}
inline PlanarLocalizer::const_iterator PlanarLocalizer::begin(void) const
    noexcept {
  return this->cbegin();
}
inline PlanarLocalizer::const_iterator PlanarLocalizer::cbegin(void) const
    noexcept {
  return reconstruction_m.cbegin();
}
inline PlanarLocalizer::iterator PlanarLocalizer::end(void) noexcept {
  return reconstruction_m.end();
}
inline PlanarLocalizer::const_iterator PlanarLocalizer::end(void) const
    noexcept {
  return this->cend();
}
inline PlanarLocalizer::const_iterator PlanarLocalizer::cend(void) const
    noexcept {
  return reconstruction_m.cend();
}

inline LargeOffsetIndex_t PlanarLocalizer::getSerializedSize(void) const {
  return reconstruction_m.getSerializedSize();
}

inline void PlanarLocalizer::serialize(ByteBuffer* a_buffer) const {
  reconstruction_m.serialize(a_buffer);
}

inline void PlanarLocalizer::unpackSerialized(ByteBuffer* a_buffer) {
  reconstruction_m.unpackSerialized(a_buffer);
}

inline PlanarLocalizer::PlanarLocalizer(const Plane& a_plane)
    : reconstruction_m(PlanarReconstructionBase::fromOnePlane(a_plane)) {}

inline PlanarLocalizer::PlanarLocalizer(const Plane& a_plane_0,
                                        const Plane& a_plane_1)
    : reconstruction_m(
          PlanarReconstructionBase::fromTwoPlanes(a_plane_0, a_plane_1)) {}

inline std::ostream& operator<<(std::ostream& out,
                                const PlanarLocalizer& a_reconstruction) {
  out << "PlanarLocalizer: \n" << a_reconstruction.reconstruction_m;
  return out;
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_LOCALIZER_TPP_
