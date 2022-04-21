// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_RECONSTRUCTION_TYPES_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_RECONSTRUCTION_TYPES_H_

#include <string>

#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_separator.h"

#include "examples/paraboloid_advector/data.h"

void getReconstruction(
    const std::string& a_reconstruction_method,
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedParaboloidLink>& a_localized_paraboloid_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::Paraboloid>* a_interface);

struct KnownCircle {
  static void getReconstruction(const Data<double>& a_liquid_volume_fraction,
                                const double a_dt, const Data<double>& a_U,
                                const Data<double>& a_V,
                                const Data<double>& a_W,
                                Data<IRL::Paraboloid>* a_interface);
};

void correctInterfacePlaneBorders(Data<IRL::Paraboloid>* a_interface);

namespace details {
inline IRL::Paraboloid fromCircle(const IRL::Pt& a_center,
                                  const double a_radius,
                                  const IRL::Normal& a_normal) {
  const double curvature = 1.0 / a_radius;
  IRL::ReferenceFrame frame;
  IRL::UnitQuaternion quat(M_PI * 0.5, IRL::Normal(0.0, 0.0, 1.0));
  frame[0] = quat * a_normal;
  frame[1] = IRL::Normal(0.0, 0.0, 1.0);
  frame[2] = a_normal;

  return IRL::Paraboloid(a_center + a_radius * a_normal, frame, 0.5 * curvature,
                         0.0);
}
}  // namespace details

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_RECONSTRUCTION_TYPES_H_
