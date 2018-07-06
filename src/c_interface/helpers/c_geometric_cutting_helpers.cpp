// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/helpers/c_geometric_cutting_helpers.h"

#include <cassert>

extern "C" {

bool c_isPtInt_PlanarSep(const double* a_pt, const c_PlanarSep* a_separator) {
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  return IRL::isPtInternal(IRL::Pt::fromRawDoublePointer(a_pt),
                           *a_separator->obj_ptr);
}

bool c_isPtInt_PlanarLoc(const double* a_pt, const c_PlanarLoc* a_localizer) {
  assert(a_localizer != nullptr);
  assert(a_localizer->obj_ptr != nullptr);
  return IRL::isPtInternal(IRL::Pt::fromRawDoublePointer(a_pt),
                           *a_localizer->obj_ptr);
}

int c_locatePt_LocSepLink(const double* a_pt,
                          const c_LocSepLink* a_locseplink) {
  assert(a_locseplink != nullptr);
  assert(a_locseplink->obj_ptr != nullptr);
  return static_cast<int>(IRL::locatePt(IRL::Pt::fromRawDoublePointer(a_pt),
                                        *a_locseplink->obj_ptr));
}

int c_locatePt_LocSepGroupLink(const double* a_pt,
                               const c_LocSepGroupLink* a_locsepgrouplink) {
  assert(a_locsepgrouplink != nullptr);
  assert(a_locsepgrouplink->obj_ptr != nullptr);
  return static_cast<int>(IRL::locatePt(IRL::Pt::fromRawDoublePointer(a_pt),
                                        *a_locsepgrouplink->obj_ptr));
}

}  // end extern C
