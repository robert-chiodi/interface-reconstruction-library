// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/parameters/constants.h"

#include <algorithm>
#include <cassert>

namespace IRL {

void setMinimumVolumeToTrack(const double a_minimum_volume_to_track) {
  global_constants::MINIMUM_VOLUME_TO_TRACK = a_minimum_volume_to_track;
}

void setMinimumSurfaceAreaToTrack(
    const double a_minimum_surface_area_to_track) {
  global_constants::MINIMUM_SURFACE_AREA_TO_TRACK =
      a_minimum_surface_area_to_track;
}

void setVolumeFractionBounds(const double a_VF_low) {
  assert(a_VF_low > 0.0 && a_VF_low < 1.0);
  global_constants::VF_LOW = a_VF_low;
  global_constants::VF_HIGH = 1.0 - a_VF_low;
  global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE =
      std::min(global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE,
               global_constants::VF_LOW);
}

void setVolumeFractionTolerance(const double a_tolerance) {
  assert(a_tolerance > 0.0 && a_tolerance < 1.0);
  global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE =
      std::min(global_constants::VF_LOW, a_tolerance);
}

namespace global_constants {

double VF_LOW = 1.0e-8;

double VF_HIGH = 1.0 - VF_LOW;

double MINIMUM_VOLUME_TO_TRACK = DBL_EPSILON;

double MINIMUM_SURFACE_AREA_TO_TRACK = std::pow(DBL_EPSILON, 2.0 / 3.0);

double TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE = 1.0e-12;

}  // namespace global_constants

}  // namespace IRL
