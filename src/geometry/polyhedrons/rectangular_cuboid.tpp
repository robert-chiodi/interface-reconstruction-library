// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_RECTANGULAR_CUBOID_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_RECTANGULAR_CUBOID_TPP_

namespace IRL {

template <class Derived, class VertexType>
double
RectangularCuboidSpecialization<Derived, VertexType>::calculateSideLength(
    const UnsignedIndex_t a_dimension) const {
  assert(a_dimension < 3);
  return ((*this)[2])[a_dimension] - ((*this)[4])[a_dimension];
}

template <class Derived, class VertexType>
inline double
RectangularCuboidSpecialization<Derived, VertexType>::calculateVolume(
    void) const {
  return this->calculateSideLength(0) * this->calculateSideLength(1) *
         this->calculateSideLength(2);
}

template <class Derived, class VertexType>
inline Pt
RectangularCuboidSpecialization<Derived, VertexType>::calculateCentroid(
    void) const {
  return 0.5 * Pt((*this)[2] + (*this)[4]);
}

template <class Derived, class VertexType>
VolumeMoments
RectangularCuboidSpecialization<Derived, VertexType>::calculateMoments(
    void) const {
  Volume volume = this->calculateVolume();
  return {volume, volume * this->calculateCentroid()};
}

template <class Derived, class VertexType>
PlanarLocalizer
RectangularCuboidSpecialization<Derived, VertexType>::getLocalizer(void) const {
  PlanarLocalizer reconstruction_to_return;
  reconstruction_to_return.setNumberOfPlanes(6);
  reconstruction_to_return[0] = Plane(Normal(-1.0, 0.0, 0.0), -(*this)[4].x());
  reconstruction_to_return[1] = Plane(Normal(1.0, 0.0, 0.0), (*this)[2].x());
  reconstruction_to_return[2] = Plane(Normal(0.0, -1.0, 0.0), -(*this)[4].y());
  reconstruction_to_return[3] = Plane(Normal(0.0, 1.0, 0.0), (*this)[2].y());
  reconstruction_to_return[4] = Plane(Normal(0.0, 0.0, -1.0), -(*this)[4].z());
  reconstruction_to_return[5] = Plane(Normal(0.0, 0.0, 1.0), (*this)[2].z());
  return reconstruction_to_return;
}

template <class VertexType>
StoredRectangularCuboid<VertexType>
StoredRectangularCuboid<VertexType>::fromBoundingPts(const Pt& a_min_point,
                                                     const Pt& a_max_point) {
  return StoredRectangularCuboid(a_min_point, a_max_point);
}

template <class VertexType>
StoredRectangularCuboid<VertexType>::StoredRectangularCuboid(
    const Pt& a_min_point, const Pt& a_max_point)
    : StoredVertexAccess<StoredRectangularCuboid<VertexType>, VertexType, 8>{
          Pt(a_max_point.x(), a_min_point.y(), a_min_point.z()),
          Pt(a_max_point.x(), a_max_point.y(), a_min_point.z()),
          Pt(a_max_point.x(), a_max_point.y(), a_max_point.z()),
          Pt(a_max_point.x(), a_min_point.y(), a_max_point.z()),
          Pt(a_min_point.x(), a_min_point.y(), a_min_point.z()),
          Pt(a_min_point.x(), a_max_point.y(), a_min_point.z()),
          Pt(a_min_point.x(), a_max_point.y(), a_max_point.z()),
          Pt(a_min_point.x(), a_min_point.y(), a_max_point.z())} {}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_RECTANGULAR_CUBOID_TPP_
