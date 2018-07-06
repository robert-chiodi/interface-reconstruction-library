// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/parameters/c_constants.h"

extern "C" {

void c_setVFBounds(const double* a_VF_low) {
  IRL::setVolumeFractionBounds(*a_VF_low);
}

void c_setVFTolerance_IterativeDistanceFinding(const double* a_tolerance) {
  IRL::setVolumeFractionTolerance(*a_tolerance);
}

void c_setMinimumVolToTrack(const double* a_minimum_volume_to_track) {
  IRL::setMinimumVolumeToTrack(*a_minimum_volume_to_track);
}

void c_setMinimumSAToTrack(const double* a_minimum_surface_area_to_track) {
  IRL::setMinimumSurfaceAreaToTrack(*a_minimum_surface_area_to_track);
}
}
