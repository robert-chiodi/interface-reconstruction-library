// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PLANE_TPP_
#define SRC_GEOMETRY_GENERAL_PLANE_TPP_

namespace IRL {

inline Plane::Plane(void) : normal_m(), distance_m{0.0} {}

inline Plane::Plane(const Normal& a_normal, const double a_distance)
    : normal_m(a_normal), distance_m(a_distance) {
  this->checkValidNormal();
}

inline Normal& Plane::normal(void) {
  this->checkValidNormal();
  return normal_m;
}

inline const Normal& Plane::normal(void) const {
  this->checkValidNormal();
  return normal_m;
}

inline double& Plane::distance(void) { return distance_m; }

inline const double& Plane::distance(void) const { return distance_m; }

inline bool Plane::operator==(const Plane& a_other_plane) const {
  return this->normal() == a_other_plane.normal() &&
         this->distance() == a_other_plane.distance();
}

inline bool Plane::operator!=(const Plane& a_other_plane) const {
  return !((*this) == a_other_plane);
}

template <class PtType>
__attribute__((const)) inline double Plane::signedDistanceToPoint(
    const PtType& a_pt) const {
  return this->normal()*a_pt.getPt() - this->distance();
}

inline LargeOffsetIndex_t Plane::getSerializedSize(void) const {
  return static_cast<LargeOffsetIndex_t>(normal_m.getSerializedSize() +
                                         sizeof(double));
}

inline void Plane::serialize(ByteBuffer* a_buffer) const {
  normal_m.serialize(a_buffer);
  a_buffer->pack(&distance_m, 1);
}

inline void Plane::unpackSerialized(ByteBuffer* a_buffer) {
  normal_m.unpackSerialized(a_buffer);
  a_buffer->unpack(&distance_m, 1);
}

inline Plane Plane::generateFlippedPlane(void) const {
  return Plane(-1.0 * normal_m, -distance_m);
}

inline void Plane::checkValidNormal(void) const {
  assert(std::fabs(normal_m * normal_m - 1.0) < 1.0e-14 ||
         normal_m * normal_m < 1.0e-14);
}

inline std::ostream& operator<<(std::ostream& out, const Plane& a_plane) {
  out << std::setprecision(15);
  out << "( " << a_plane.normal()[0];
  out << ", " << a_plane.normal()[1];
  out << ", " << a_plane.normal()[2];
  out << ", " << a_plane.distance();
  out << " ) \n";
  return out;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_PLANE_TPP_
