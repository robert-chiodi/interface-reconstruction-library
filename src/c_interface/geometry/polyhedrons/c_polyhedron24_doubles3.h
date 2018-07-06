// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_POLYHEDRON24_DOUBLES3_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_POLYHEDRON24_DOUBLES3_H_

#include "src/geometry/general/new_pt_calculation_functors.h"
#include "src/geometry/general/pt_with_data.h"
#include "src/geometry/polyhedrons/polyhedron_24.h"

extern "C" {

struct c_Poly24_d3 {
  IRL::StoredPolyhedron24<
      IRL::PtWithDoublesStatelessFunctor<IRL::LinearInterpolation_Functor, 3>>*
      obj_ptr = nullptr;
};

void c_Poly24_d3_new(c_Poly24_d3* a_ptr);

void c_Poly24_d3_delete(c_Poly24_d3* a_ptr);

void c_Poly24_d3_construct(c_Poly24_d3* a_ptr, const double* a_dodecahedron,
                           const double* a_attached_data);

void c_Poly24_d3_adjustCapToMatchVolume(c_Poly24_d3* a_ptr,
                                        const double* a_correct_signed_volume);

void c_Poly24_d3_getBoundingPts(c_Poly24_d3* a_ptr, double* a_lower_pt,
                                double* a_upper_pt);

void c_Poly24_d3_getPt(c_Poly24_d3* a_ptr, const int* a_index, double* a_pt);

void c_Poly24_d3_setPt(c_Poly24_d3* a_ptr, const int* a_index,
                       const double* a_pt);

void c_Poly24_d3_getData(c_Poly24_d3* a_ptr, const int* a_index,
                         double* a_data);

void c_Poly24_d3_setData(c_Poly24_d3* a_ptr, const int* a_index,
                         const double* a_data);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_POLYHEDRON24_DOUBLES3_H_
