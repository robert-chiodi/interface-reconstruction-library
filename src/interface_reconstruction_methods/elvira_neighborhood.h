// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_ELVIRA_NEIGHBORHOOD_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_ELVIRA_NEIGHBORHOOD_H_

#include <float.h>

#include <cassert>
#include <string>

#include "src/generic_cutting/cut_polygon.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/moments/cell_collection.h"
#include "src/moments/cell_grouped_moments.h"
#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief Below plane volume fraction information in a stencil, to be used
/// for ELVIRA. Up to 27 cells to cover 3x3x3 stencil in 3D. Can partially
/// fill to 9 for 2D.
///
/// This class is used to store neighboring volume fraction information and
/// RectangularCuboid cells to be used with the ELVIRA. There is an assumed
/// order of the cells that will be added, with x being the fastest dimension, y
/// being second, and z being slowest. i.e. for the loop order:
/// <br></br> loop z
///   loop y
///     loop x
///       (*this)[] = (...)
class ELVIRANeighborhood {
  using CGD = CellGroupedMoments<RectangularCuboid, double>;

 public:
  /// \brief Default constructor.
  ELVIRANeighborhood(void) = default;

  /// \brief Construct a CellGroupedMoments and add it
  /// to the collection for index i,j,k.
  void setMember(const RectangularCuboid* a_rectangular_cuboid,
                 const double* a_liquid_volume_fraction, const int i,
                 const int j, const int k = -1);

  /// \brief Return the cell stored at the index i,j,k
  const RectangularCuboid& getCell(const int i, const int j,
                                   const int k = -1) const;

  /// \brief Return moments stored at the index i,j,k
  double getStoredMoments(const int i, const int j, const int k = -1) const;

  /// \brief Set size of the neighborhood.
  void resize(const UnsignedIndex_t a_size);

  /// \brief Default destructor.
  ~ELVIRANeighborhood(void) = default;

 private:
  /// \brief Calculate linear index from i,j,k.
  UnsignedIndex_t calculateLinearIndex(const int i, const int j,
                                       const int k) const;

  CellCollection<CGD>
      collection_m;  ///< \brief Collection that holds correct moments.
};

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_ELVIRA_NEIGHBORHOOD_H_
