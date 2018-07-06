// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_NEW_PT_CALCULATION_FUNCTORS_H_
#define SRC_GEOMETRY_GENERAL_NEW_PT_CALCULATION_FUNCTORS_H_

namespace IRL {

class LinearInterpolation_Functor {
 public:
  template <class PtWithDataType>
  PtWithDataType operator()(const PtWithDataType& a_pt_0, const double a_dist_0,
                            const PtWithDataType& a_pt_1,
                            const double a_dist_1);
};

}  // namespace IRL

#include "src/geometry/general/new_pt_calculation_functors.tpp"

#endif  // SRC_GEOMETRY_GENERAL_NEW_PT_CALCULATION_FUNCTORS_H_
