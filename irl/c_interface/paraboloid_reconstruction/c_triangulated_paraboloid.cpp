// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/paraboloid_reconstruction/c_triangulated_paraboloid.h"

#include <iostream>

extern "C" {

void c_TriangulatedParaboloid_new(c_TriangulatedParaboloid* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::TriangulatedSurfaceOutput;
}

void c_TriangulatedParaboloid_newFromObjectAllocationServer(
    c_TriangulatedParaboloid* a_self,
    c_ObjServer_TriangulatedParaboloid* a_object_allocation_server) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
}

void c_TriangulatedParaboloid_delete(c_TriangulatedParaboloid* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_TriangulatedParaboloid_clear(c_TriangulatedParaboloid* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->clearAll();
}

int c_TriangulatedParaboloid_getNumberOfTriangles(
    const c_TriangulatedParaboloid* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  // std::cout << "Surface has " << a_self->obj_ptr->nTriangles() << "
  // triangles"
  //           << std::endl;
  return static_cast<int>(a_self->obj_ptr->nTriangles());
}

void c_TriangulatedParaboloid_getPt(c_TriangulatedParaboloid* a_self,
                                    const int* a_tri_index,
                                    const int* a_vert_index, double* a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_tri_index >= 0 &&
         *a_tri_index < static_cast<int>(a_self->obj_ptr->nTriangles()));
  assert(*a_vert_index >= 0 && *a_vert_index < 3);
  const auto& triangle =
      a_self->obj_ptr
          ->getTriangleList()[static_cast<IRL::UnsignedIndex_t>(*a_tri_index)];
  // std::cout << "Pt[" << *a_tri_index << "][" << *a_vert_index << "] = ";
  // std::cout << triangle[static_cast<IRL::UnsignedIndex_t>(*a_vert_index)][0]
  //           << ", ";
  // std::cout << triangle[static_cast<IRL::UnsignedIndex_t>(*a_vert_index)][1]
  //           << ", ";
  // std::cout << triangle[static_cast<IRL::UnsignedIndex_t>(*a_vert_index)][2]
  //           << std::endl;
  a_pt[0] = triangle[static_cast<IRL::UnsignedIndex_t>(*a_vert_index)][0];
  a_pt[1] = triangle[static_cast<IRL::UnsignedIndex_t>(*a_vert_index)][1];
  a_pt[2] = triangle[static_cast<IRL::UnsignedIndex_t>(*a_vert_index)][2];
}

}  // end extern C
