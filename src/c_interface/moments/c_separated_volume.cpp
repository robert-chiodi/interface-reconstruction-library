// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_separated_volume.h"

extern "C" {

void c_SepVol_new(c_SepVol* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::SeparatedMoments<IRL::Volume>;
}

void c_SepVol_delete(c_SepVol* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_SepVol_construct(c_SepVol* a_self, const double* a_moments_list) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr) =
      IRL::SeparatedMoments<IRL::Volume>::fromRawDoublePointer(
          a_moments_list);
}

double c_SepVol_getVolume(const c_SepVol* a_self, const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < 2);
  return (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)];
}


}
