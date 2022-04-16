// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_REFERENCE_FRAME_H_
#define IRL_GEOMETRY_GENERAL_REFERENCE_FRAME_H_

#include <array>

#include "irl/geometry/general/normal.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief A reference frame with three normals (1 for each direction in 3D
/// space).
class ReferenceFrame {
  using iterator = typename std::array<Normal, 3>::iterator;
  using const_iterator = typename std::array<Normal, 3>::const_iterator;

 public:
  /// \brief Default constructor.
  ReferenceFrame(void) = default;

  /// \brief Construct given 3 normals
  ReferenceFrame(const Normal& a_axis_0, const Normal& a_axis_1,
                 const Normal& a_axis_2);

  /// \brief Overload `operator[]` for access.
  Normal& operator[](const UnsignedIndex_t a_axis);

  /// \brief Const version of verload `operator[]` for access.
  const Normal& operator[](const UnsignedIndex_t a_axis) const;

  bool isOrthonormalBasis(void) const;

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~ReferenceFrame(void) = default;

 private:
  /// \brief Three orthonormal vectors making up the reference frame.
  std::array<Normal, 3> axis_m;
};

}  // namespace IRL

#include "irl/geometry/general/reference_frame.tpp"

#endif  // IRL_GEOMETRY_GENERAL_REFERENCE_FRAME_H_
