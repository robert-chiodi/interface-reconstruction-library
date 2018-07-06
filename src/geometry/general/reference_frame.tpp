// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_REFERENCE_FRAME_TPP_
#define SRC_GEOMETRY_GENERAL_REFERENCE_FRAME_TPP_

namespace IRL {

inline ReferenceFrame::ReferenceFrame(const Normal& a_axis_0,
                                      const Normal& a_axis_1,
                                      const Normal& a_axis_2)
    : axis_m{a_axis_0, a_axis_1, a_axis_2} {}

inline Normal& ReferenceFrame::operator[](const UnsignedIndex_t a_axis) {
  assert(a_axis < 3);
  return axis_m[a_axis];
}

inline const Normal& ReferenceFrame::operator[](
    const UnsignedIndex_t a_axis) const {
  assert(a_axis < 3);
  return axis_m[a_axis];
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_REFERENCE_FRAME_TPP_
