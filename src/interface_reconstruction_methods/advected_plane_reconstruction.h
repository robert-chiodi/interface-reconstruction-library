// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_ADVECTED_PLANE_RECONSTRUCTION_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_ADVECTED_PLANE_RECONSTRUCTION_H_

#include <string>

#include "src/distributions/k_means.h"
#include "src/distributions/partition_by_normal_vector.h"
#include "src/generic_cutting/cut_polygon.h"
#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/plane.h"
#include "src/helpers/mymath.h"
#include "src/interface_reconstruction_methods/plane_distance.h"
#include "src/interface_reconstruction_methods/r2p_optimization.h"
#include "src/interface_reconstruction_methods/reconstruction_cleaning.h"
#include "src/moments/separated_volume_moments.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

class AdvectedPlaneReconstruction {
  friend class AdvectedPlaneReconstructionDebug;

 public:
  /// \brief Solve for the AdvectedPlaneReconstruction, returning a
  /// PlanarSeparator.
  template <class MomentsContainerType, class CellType>
  static PlanarSeparator solve(
      const MomentsContainerType& a_volume_moments_list,
      const R2PNeighborhood<CellType>& a_neighborhood,
      const double threshold = 0.9);

 private:
  template <class ContainedType, class CellType>
  static PlanarSeparator findBestPermutation(
      ContainedType* a_moments, const R2PNeighborhood<CellType>& a_neighborhood,
      const double a_target_volume_fraction);

  template <class CellType>
  static PlanarSeparator constructSeparatorAttempt(
      const CellType& a_cell, const double a_target_volume_fraction,
      const Normal& a_normal_0, const Pt& a_pt_0, const Normal& a_normal_1,
      const Pt& a_pt_1, const double a_flip_cut);

  template <class CellType>
  static void checkIfBest(const R2PNeighborhood<CellType>& a_neighborhood,
                          const PlanarSeparator& a_attempt_separator,
                          PlanarSeparator* a_current_best_separator,
                          double* a_current_minimum_error);
};

class AdvectedPlaneReconstructionDebug {
 public:
  template <class MomentsContainerType, class CellType>
  static PlanarSeparator solve(
      const MomentsContainerType& a_volume_moments_list,
      const R2PNeighborhood<CellType>& a_neighborhood,
      const double threshold = 0.9);

 private:
  template <class ContainedType, class CellType>
  static PlanarSeparator findBestPermutation(
      ContainedType* a_moments, const R2PNeighborhood<CellType>& a_neighborhood,
      const double a_target_volume_fraction);

  template <class CellType>
  static void checkIfBest(const R2PNeighborhood<CellType>& a_neighborhood,
                          const PlanarSeparator& a_attempt_separator,
                          PlanarSeparator* a_current_best_separator,
                          double* a_current_minimum_error);

  template <class CellType>
  static void writeOutPlane(const R2PNeighborhood<CellType>& a_neighborhood,
                            const PlanarSeparator& a_reconstruction,
                            const std::string& a_prefix,
                            const std::size_t a_iteration_number);

  template <class CellType>
  static void writeOutCentroids(
      const R2PNeighborhood<CellType>& a_neighborhood);
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/advected_plane_reconstruction.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_ADVECTED_PLANE_RECONSTRUCTION_H_
