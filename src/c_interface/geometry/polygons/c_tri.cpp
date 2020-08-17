// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polygons/c_tri.h"

#include <cassert>

extern "C" {

void c_Tri_new(c_Tri* a_ptr) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr == nullptr);
  a_ptr->obj_ptr = new IRL::Tri;
}

void c_Tri_delete(c_Tri* a_ptr) {
  delete a_ptr->obj_ptr;
  a_ptr->obj_ptr = nullptr;
}

void c_Tri_construct(c_Tri* a_ptr, const double* a_pts) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  (*a_ptr->obj_ptr) = IRL::Tri::fromRawDoublePointer(3, a_pts);
}

void c_Tri_getVertices(c_Tri* a_ptr, double* a_pts) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    a_pts[3 * n + 0] = (*a_ptr->obj_ptr)[n].x();
    a_pts[3 * n + 1] = (*a_ptr->obj_ptr)[n].y();
    a_pts[3 * n + 2] = (*a_ptr->obj_ptr)[n].z();
  }
}

double c_Tri_calculateVolume(c_Tri* a_ptr) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  return a_ptr->obj_ptr->calculateVolume();
}

void c_Tri_calculateCentroid(c_Tri* a_ptr, double* a_centroid) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  IRL::Pt centroid = a_ptr->obj_ptr->calculateCentroid();
  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    a_centroid[n] = centroid[n];
  }
}

void c_Tri_calculateNormal(c_Tri* a_ptr, double* a_normal) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  IRL::Normal face_normal(a_ptr->obj_ptr->calculateNormal());
  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    a_normal[n] = face_normal[n];
  }
}

void c_Tri_getLocalizer(const c_Tri* a_ptr, c_PlanarLoc* a_planar_localizer) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  assert(a_planar_localizer != nullptr);
  assert(a_planar_localizer->obj_ptr != nullptr);
  *a_planar_localizer->obj_ptr = a_ptr->obj_ptr->getLocalizer();
}

void c_Tri_reversePtOrdering(c_Tri* a_ptr) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  a_ptr->obj_ptr->reversePtOrdering();
}

void c_Tri_getBoundingPts(c_Tri* a_ptr, double* __restrict__ a_lower_pt,
                          double* __restrict__ a_upper_pt) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  IRL::Pt lower_pt = a_ptr->obj_ptr->getLowerLimits();
  IRL::Pt upper_pt = a_ptr->obj_ptr->getUpperLimits();
  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    a_lower_pt[n] = lower_pt[n];
    a_upper_pt[n] = upper_pt[n];
  }
}

double c_Tri_calculateSign(const c_Tri* a_ptr) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  return a_ptr->obj_ptr->calculateSign();
}

void c_Tri_setPlaneOfExistence(c_Tri* a_ptr, const double* a_plane) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  a_ptr->obj_ptr->setPlaneOfExistence(
      IRL::Plane(IRL::Normal(a_plane[0], a_plane[1], a_plane[2]), a_plane[3]));
}

void c_Tri_calculateAndSetPlaneOfExistence(c_Tri* a_ptr) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  a_ptr->obj_ptr->calculateAndSetPlaneOfExistence();
}

void c_Tri_getPlaneOfExistence(const c_Tri* a_ptr, double* a_plane) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  const IRL::Plane& plane = a_ptr->obj_ptr->getPlaneOfExistence();
  a_plane[0] = plane.normal()[0];
  a_plane[1] = plane.normal()[1];
  a_plane[2] = plane.normal()[2];
  a_plane[3] = plane.distance();
}
}
