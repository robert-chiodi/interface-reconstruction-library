// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_CELL_GROUPED_MOMENTS_H_
#define SRC_MOMENTS_CELL_GROUPED_MOMENTS_H_

#include "src/generic_cutting/generic_cutting.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

/// \brief A class that couples together a cell and another class.
template <class CellType, class ContainedMomentsType>
class CellGroupedMoments {
 public:
  using cell_type = CellType;
  using contained_type = ContainedMomentsType;

  /// \brief Default constructor.
  CellGroupedMoments(void);

  /// \brief Construct given pointers to `a_cell` and `a_contained_type`.
  CellGroupedMoments(const CellType* a_cell,
                     const ContainedMomentsType* a_contained_object);

  /// \brief Return const reference to the pointed to cell.
  const CellType& getCell(void) const;

  /// \brief Return const reference to the pointed to cell.
  const ContainedMomentsType& getStoredMoments(void) const;

  /// \brief Return the normalized ContainedType for the cell given
  /// `a_separator`.
  ContainedMomentsType calculateNormalizedVolumeMoments(
      const PlanarSeparator& a_separator) const;

  /// \brief Default destructor.
  ~CellGroupedMoments(void) = default;

 private:
  const CellType* ptr_to_cell_m;
  const ContainedMomentsType* ptr_to_contained_m;
};

}  // namespace IRL

#include "src/moments/cell_grouped_moments.tpp"

#endif  // SRC_MOMENTS_CELL_GROUPED_MOMENTS_H_
