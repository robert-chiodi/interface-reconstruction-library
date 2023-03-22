// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_PLANE_TPP_
#define IRL_GEOMETRY_GENERAL_PLANE_TPP_

namespace IRL {

template <class ScalarType>
inline PlaneBase<ScalarType>::PlaneBase(void)
    : normal_m(), distance_m(static_cast<ScalarType>(0.0)) {}

template <class ScalarType>
inline PlaneBase<ScalarType>::PlaneBase(const NormalBase<ScalarType>& a_normal,
                                        const ScalarType a_distance)
    : normal_m(a_normal), distance_m(a_distance) {
  this->checkValidNormal();
}

template <class ScalarType>
inline NormalBase<ScalarType>& PlaneBase<ScalarType>::normal(void) {
  this->checkValidNormal();
  return normal_m;
}

template <class ScalarType>
inline const NormalBase<ScalarType>& PlaneBase<ScalarType>::normal(void) const {
  this->checkValidNormal();
  return normal_m;
}

template <class ScalarType>
inline ScalarType& PlaneBase<ScalarType>::distance(void) {
  return distance_m;
}

template <class ScalarType>
inline const ScalarType& PlaneBase<ScalarType>::distance(void) const {
  return distance_m;
}

template <class ScalarType>
inline bool PlaneBase<ScalarType>::operator==(
    const PlaneBase<ScalarType>& a_other_plane) const {
  return this->normal() == a_other_plane.normal() &&
         this->distance() == a_other_plane.distance();
}

template <class ScalarType>
inline bool PlaneBase<ScalarType>::operator!=(
    const PlaneBase<ScalarType>& a_other_plane) const {
  return !((*this) == a_other_plane);
}

template <class ScalarType>
template <class PtType>
__attribute__((const)) inline ScalarType
PlaneBase<ScalarType>::signedDistanceToPoint(const PtType& a_pt) const {
  return this->normal() * a_pt.getPt() - this->distance();
}

template <class ScalarType>
inline LargeOffsetIndex_t PlaneBase<ScalarType>::getSerializedSize(void) const {
  return static_cast<LargeOffsetIndex_t>(normal_m.getSerializedSize() +
                                         sizeof(ScalarType));
}

template <class ScalarType>
inline void PlaneBase<ScalarType>::serialize(ByteBuffer* a_buffer) const {
  normal_m.serialize(a_buffer);
  a_buffer->pack(&distance_m, 1);
}

template <class ScalarType>
inline void PlaneBase<ScalarType>::unpackSerialized(ByteBuffer* a_buffer) {
  normal_m.unpackSerialized(a_buffer);
  a_buffer->unpack(&distance_m, 1);
}

template <class ScalarType>
inline PlaneBase<ScalarType> PlaneBase<ScalarType>::generateFlippedPlane(
    void) const {
  return PlaneBase<ScalarType>(-static_cast<ScalarType>(1) * normal_m,
                               -distance_m);
}

template <class ScalarType>
inline void PlaneBase<ScalarType>::checkValidNormal(void) const {
  assert(fabs(normal_m * normal_m - static_cast<ScalarType>(1)) <
             static_cast<ScalarType>(1.0e-14) ||
         normal_m * normal_m < static_cast<ScalarType>(1.0e-14));
}

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const PlaneBase<ScalarType>& a_plane) {
  out << std::setprecision(15);
  out << "( " << a_plane.normal()[0];
  out << ", " << a_plane.normal()[1];
  out << ", " << a_plane.normal()[2];
  out << ", " << a_plane.distance();
  out << " ) \n";
  return out;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_PLANE_TPP_
