// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_NEIGHBORHOOD_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_NEIGHBORHOOD_H_

#include "src/moments/cell_collection.h"
#include "src/moments/cell_grouped_moments.h"

namespace IRL {

/// \brief Neighborhood storage used in the R2P
/// optimization routines.
template <class CellType>
class R2PNeighborhood {
  using CGD = CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>;
  using iterator = typename CellCollection<CGD>::iterator;
  using const_iterator = typename CellCollection<CGD>::const_iterator;

 public:
  using cell_type = CellType;

  /// \brief Default constructor.
  R2PNeighborhood(void);

  /// \brief Construct CellGroupedMoments and add to end of collection.
  void addMember(const CellType* a_cell,
                 const SeparatedMoments<VolumeMoments>* a_volume_moments);
  /// \brief Reset neighborhood size to 0.
  void emptyNeighborhood(void);

  /// \brief Construct CellGroupedMoments and place into collection.
  void setMember(const UnsignedIndex_t a_index, const CellType* a_cell,
                 const SeparatedMoments<VolumeMoments>* a_volume_moments);

  /// \brief Set the index for the center cell in the collection.
  void setCenterOfStencil(const UnsignedIndex_t a_index);

  /// \brief Set the index for the center cell in the collection.
  void setSurfaceArea(const double a_surface_area);

  /// \brief Return the center cell.
  const CellType& getCenterCell(void) const;

  /// \brief Return the center cell moments
  const SeparatedMoments<VolumeMoments>& getCenterCellStoredMoments(void) const;

  /// \brief Return the cell stored at the index
  const typename CGD::cell_type& getCell(const UnsignedIndex_t a_index) const;

  /// \brief Return moments stored at the index
  const SeparatedMoments<VolumeMoments>& getStoredMoments(
      const UnsignedIndex_t a_index) const;

  /// \brief Set size of the neighborhood.
  void resize(const UnsignedIndex_t a_size);

  /// \brief Get size of the collection.
  UnsignedIndex_t size(void) const;

  /// \brief Get surface area for the center cell
  double getSurfaceArea(void) const;

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~R2PNeighborhood(void) = default;

 private:
  /// \brief Make sure index is not larger than current collection size.
  void checkIndex(UnsignedIndex_t a_index) const;
  void checkCenterStencilSet(void) const;

  /// \brief Collection of cells and cell moments.
  CellCollection<CGD> collection_m;
  /// brief Interface surface area for center cell in neighborhood.
  double center_cell_surface_area_m;
  /// \brief Center stencil cell index in the list of added cells.
  UnsignedIndex_t center_cell_index_m;
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/r2p_neighborhood.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_NEIGHBORHOOD_H_
