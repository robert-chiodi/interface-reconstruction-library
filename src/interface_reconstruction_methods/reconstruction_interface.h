// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_INTERFACE_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_INTERFACE_H_

#include "src/interface_reconstruction_methods/advected_plane_reconstruction.h"
#include "src/interface_reconstruction_methods/elvira.h"
#include "src/interface_reconstruction_methods/lvira_neighborhood.h"
#include "src/interface_reconstruction_methods/lvira_optimization.h"
#include "src/interface_reconstruction_methods/mof.h"
#include "src/interface_reconstruction_methods/r2p_optimization.h"
#include "src/interface_reconstruction_methods/reconstruction_cleaning.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

/// \brief Perform R2P reconstruction for a 2D problem in the x-y plane.
template <class CellType>
inline PlanarSeparator reconstructionWithR2P2D(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform R2P reconstruction for a 3D problem.
template <class CellType>
inline PlanarSeparator reconstructionWithR2P3D(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform ELVIRA Reconstruction for 2D.
inline PlanarSeparator reconstructionWithELVIRA2D(
    const ELVIRANeighborhood& a_neighborhood_geometry);

/// \brief Perform ELVIRA Reconstruction for 3D.
inline PlanarSeparator reconstructionWithELVIRA3D(
    const ELVIRANeighborhood& a_neighborhood_geometry);

/// \brief Perform LVIRA Reconstruction for 2D.
template <class CellType>
inline PlanarSeparator reconstructionWithELVIRA2D(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform LVIRA Reconstruction for 3D.
template <class CellType>
inline PlanarSeparator reconstructionWithELVIRA3D(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform MOF Reconstruction for 2D with optional weights.
/// Defaults to even weighting.
template <class CellType>
PlanarSeparator reconstructionWithMOF2D(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight = 0.5, const double a_external_weight = 0.5);

/// \brief Perform MOF Reconstruction for 3D with optional weights.
/// Defaults to even weighting.
template <class CellType>
PlanarSeparator reconstructionWithMOF3D(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight = 0.5, const double a_external_weight = 0.5);

/// \brief Reconstruct a PlanarSeparator from a collection of
/// VolumeMomentsAndNormal objects, with supplying the
/// threshold for when one or two planes are used.
template <class MomentsContainerType, class CellType>
PlanarSeparator reconstructionWithAdvectedNormals(
    const MomentsContainerType& a_volume_moments_list,
    const R2PNeighborhood<CellType>& a_neighborhood,
    const double a_two_plane_threshold = 0.95);

//******************************************************************* //
//     Debug versions below this, which export solution
//      process to screen at end.
//******************************************************************* //
/// \brief Perform R2P reconstruction for a 2D problem in the x-y plane
/// and prints out all accepted reconstructions on path to best solution.
template <class CellType>
inline PlanarSeparator reconstructionWithR2P2DDebug(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform R2P reconstruction for a 3D problem
/// and prints out all accepted reconstructions on path to best solution.
template <class CellType>
inline PlanarSeparator reconstructionWithR2P3DDebug(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform ELVIRA Reconstruction for 2D and
/// prints out each of attempted reconstructions.
inline PlanarSeparator reconstructionWithELVIRA2DDebug(
    const ELVIRANeighborhood& a_neighborhood_geometry);

/// \brief Perform ELVIRA Reconstruction for 3D and
/// prints out each of attempted reconstructions.
inline PlanarSeparator reconstructionWithELVIRA3DDebug(
    const ELVIRANeighborhood& a_neighborhood_geometry);

/// \brief Perform LVIRA reconstruction for a 2D problem in the x-y plane
/// and prints out all accepted reconstructions on path to best solution.
template <class CellType>
inline PlanarSeparator reconstructionWithLVIRA2DDebug(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform LVIRA reconstruction for a 3D problem
/// and prints out all accepted reconstructions on path to best solution.
template <class CellType>
inline PlanarSeparator reconstructionWithLVIRA3DDebug(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction);

/// \brief Perform MOF Reconstruction for a 2D problem
/// and prints out all accepted reconstructions on path to best solution.
template <class CellType>
PlanarSeparator reconstructionWithMOF2DDebug(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight = 0.5, const double a_external_weight = 0.5);

/// \brief Perform MOF Reconstruction for a 3D problem
/// and prints out all accepted reconstructions on path to best solution.
template <class CellType>
PlanarSeparator reconstructionWithMOF3DDebug(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight = 0.5, const double a_external_weight = 0.5);

/// \brief Reconstruct a PlanarSeparator from a collection of
/// VolumeMomentsAndNormal objects, with supplying the
/// threshold for when one or two planes are used.
template <class MomentsContainerType, class CellType>
PlanarSeparator reconstructionWithAdvectedNormalsDebug(
    const MomentsContainerType& a_volume_moments_list,
    const R2PNeighborhood<CellType>& a_neighborhood,
    const double a_two_plane_threshold = 0.95);

}  // namespace IRL

#include "src/interface_reconstruction_methods/reconstruction_interface.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_INTERFACE_H_
