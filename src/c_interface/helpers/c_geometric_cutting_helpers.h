// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_HELPERS_C_GEOMETRIC_CUTTING_HELPERS_H_
#define SRC_C_INTERFACE_HELPERS_C_GEOMETRIC_CUTTING_HELPERS_H_

#include "src/c_interface/planar_reconstruction/c_localizers.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/c_interface/planar_reconstruction/c_localized_separator_link.h"
#include "src/c_interface/planar_reconstruction/c_localized_separator_group_link.h"
#include "src/helpers/geometric_cutting_helpers.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/planar_separator.h"

extern "C" {

bool c_isPtInt_PlanarSep(const double* a_pt, const c_PlanarSep* a_separator);

bool c_isPtInt_PlanarLoc(const double* a_pt, const c_PlanarLoc* a_localizer);

int c_locatePt_LocSepLink(const double* a_pt, const c_LocSepLink* a_locseplink);

int c_locatePt_LocSepGroupLink(const double* a_pt, const c_LocSepGroupLink* a_locsepgrouplink);

}  // end extern C

#endif  // SRC_C_INTERFACE_HELPERS_C_GEOMETRIC_CUTTING_HELPERS_H_
