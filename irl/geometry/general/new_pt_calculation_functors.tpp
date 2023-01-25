// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_NEW_PT_CALCULATION_FUNCTORS_TPP_
#define IRL_GEOMETRY_GENERAL_NEW_PT_CALCULATION_FUNCTORS_TPP_

#include "irl/helpers/helper.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class PtWithDataType>
PtWithDataType LinearInterpolation_Functor::operator()(
    const PtWithDataType& a_pt_0,
    const typename PtWithDataType::value_type a_dist_0,
    const PtWithDataType& a_pt_1,
    const typename PtWithDataType::value_type a_dist_1) {
  using ScalarType = typename PtWithDataType::value_type;
  const ScalarType mu_1 = -a_dist_0 / safelyTiny(a_dist_1 - a_dist_0);
  const ScalarType mu_0 = 1.0 - mu_1;
  const auto& base_pt_0 = a_pt_0.getPt();
  const auto& base_pt_1 = a_pt_1.getPt();
  const auto& data_0 = a_pt_0.getData();
  const auto& data_1 = a_pt_1.getData();

  PtWithDataType pt_to_return;
  pt_to_return.getPt() = mu_0 * base_pt_0 + mu_1 * base_pt_1;
  auto& new_data = pt_to_return.getData();
  for (UnsignedIndex_t n = 0; n < data_0.size(); ++n) {
    new_data[n] = mu_0 * data_0[n] + mu_1 * data_1[n];
  }
  return pt_to_return;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_NEW_PT_CALCULATION_FUNCTORS_TPP_
