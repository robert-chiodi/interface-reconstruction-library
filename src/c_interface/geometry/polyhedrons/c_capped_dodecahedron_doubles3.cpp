// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/c_capped_dodecahedron_doubles3.h"

#include <cassert>

extern "C" {

void c_CapDod_d3_new(c_CapDod_d3* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::StoredCappedDodecahedron<
      IRL::PtWithDoublesStatelessFunctor<IRL::LinearInterpolation_Functor, 3>>;
}

void c_CapDod_d3_delete(c_CapDod_d3* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_CapDod_d3_construct(c_CapDod_d3* a_self,
                           const double* a_capped_dodecahedron,
                           const double* a_attached_data) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  for (IRL::UnsignedIndex_t n = 0; n < a_self->obj_ptr->getNumberOfVertices();
       ++n) {
    (*a_self->obj_ptr)[n] =
        IRL::PtWithDoublesStatelessFunctor<IRL::LinearInterpolation_Functor, 3>(
            IRL::Pt::fromRawDoublePointer(&a_capped_dodecahedron[3 * n]),
            std::array<double, 3>{a_attached_data[3 * n],
                                  a_attached_data[3 * n + 1],
                                  a_attached_data[3 * n + 2]});
  }
}

void c_CapDod_d3_adjustCapToMatchVolume(c_CapDod_d3* a_self,
                                        const double* a_correct_signed_volume) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->adjustCapToMatchVolume(*a_correct_signed_volume);
}

void c_CapDod_d3_getBoundingPts(const c_CapDod_d3* a_self, double* a_lower_pt,
                                double* a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  const IRL::Pt& lower_pt = a_self->obj_ptr->getLowerLimits();
  const IRL::Pt& upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}

void c_CapDod_d3_getPt(const c_CapDod_d3* a_self, const int* a_index,
                       double* a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  a_pt[0] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                .getPt()
                .x();
  a_pt[1] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                .getPt()
                .y();
  a_pt[2] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                .getPt()
                .z();
}

void c_CapDod_d3_setPt(c_CapDod_d3* a_self, const int* a_index,
                       const double* a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)][0] = a_pt[0];
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)][1] = a_pt[1];
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)][2] = a_pt[2];
}

void c_CapDod_d3_getData(c_CapDod_d3* a_self, const int* a_index,
                         double* a_data) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  a_data[0] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                  .getData()[0];
  a_data[1] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                  .getData()[1];
  a_data[2] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                  .getData()[2];
}

void c_CapDod_d3_setData(c_CapDod_d3* a_self, const int* a_index,
                         const double* a_data) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].getData()[0] =
      a_data[0];
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].getData()[1] =
      a_data[1];
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].getData()[2] =
      a_data[2];
}
}
