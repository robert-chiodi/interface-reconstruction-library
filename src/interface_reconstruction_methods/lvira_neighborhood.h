// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_NEIGHBORHOOD_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_NEIGHBORHOOD_H_

#include "src/moments/cell_collection.h"
#include "src/moments/cell_grouped_moments.h"

namespace IRL {

/// \brief Neighborhood storage used in the LVIRA
/// optimization routines. This stores the CellGroupedMoments
/// of the cell and the volume fraction.
template <class CellType>
class LVIRANeighborhood {
  using CGD = CellGroupedMoments<CellType, double>;
  using iterator = typename CellCollection<CGD>::iterator;
  using const_iterator = typename CellCollection<CGD>::const_iterator;

 public:
  using cell_type = CellType;

  /// \brief Default constructor.
  LVIRANeighborhood(void);

  /// \brief Construct CellGroupedMoments and add to end of collection.
  void addMember(const CellType* a_cell, const double* a_volume_fraction);

  /// \brief Construct CellGroupedMoments and place into collection.
  void setMember(const UnsignedIndex_t a_index, const CellType* a_cell,
                 const double* a_volume_fraction);

  /// \brief Reset neighborhood size to 0.
  void emptyNeighborhood(void);

  /// \brief Set size of the neighborhood.
  void resize(const UnsignedIndex_t a_size);

  /// \brief Set the index for the center cell in the collection.
  void setCenterOfStencil(const UnsignedIndex_t a_index);

  /// \brief Return the index for the center stencil
  UnsignedIndex_t getCenterOfStencilIndex(void) const;

  /// \brief Return the center cell.
  const CellType& getCenterCell(void) const;

  /// \brief Return the center cell moments
  const double& getCenterCellStoredMoments(void) const;

  /// \brief Return the cell stored at the index
  const typename CGD::cell_type& getCell(const UnsignedIndex_t a_index) const;

  /// \brief Return moments stored at the index
  const double& getStoredMoments(const UnsignedIndex_t a_index) const;

  /// \brief Get size of the collection.
  UnsignedIndex_t size(void) const;

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~LVIRANeighborhood(void) = default;

 private:
  /// \brief Make sure index is not larger than current collection size.
  void checkIndex(UnsignedIndex_t a_index) const;
  void checkCenterStencilSet(void) const;

  /// \brief Collection of cells and cell moments.
  CellCollection<CGD> collection_m;
  /// \brief Center stencil cell index in the list of added cells.
  UnsignedIndex_t center_cell_index_m;
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/lvira_neighborhood.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_NEIGHBORHOOD_H_
