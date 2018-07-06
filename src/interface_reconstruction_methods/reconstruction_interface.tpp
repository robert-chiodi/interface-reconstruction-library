// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_INTERFACE_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_INTERFACE_TPP_

namespace IRL {

template <class CellType>
PlanarSeparator reconstructionWithR2P2D(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  cleanReconstruction(
      a_neighborhood_geometry.getCenterCell(),
      (a_neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          a_neighborhood_geometry.getCenterCell().calculateVolume(),
      &a_initial_reconstruction);
  if (a_initial_reconstruction.getNumberOfPlanes() == 1) {
    R2P_2D1P<CellType> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  } else {
    R2P_2D2P<CellType> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  }
}

template <class CellType>
PlanarSeparator reconstructionWithR2P3D(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  cleanReconstruction(
      a_neighborhood_geometry.getCenterCell(),
      (a_neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          a_neighborhood_geometry.getCenterCell().calculateVolume(),
      &a_initial_reconstruction);
  if (a_initial_reconstruction.getNumberOfPlanes() == 1) {
    R2P_3D1P<CellType> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  } else {
    R2P_3D2P<CellType> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  }
}

PlanarSeparator reconstructionWithELVIRA2D(
    const ELVIRANeighborhood& a_neighborhood_geometry) {
  ELVIRA_2D elvira_system;
  return elvira_system.solve(&a_neighborhood_geometry);
}

PlanarSeparator reconstructionWithELVIRA3D(
    const ELVIRANeighborhood& a_neighborhood_geometry) {
  ELVIRA_3D elvira_system;
  return elvira_system.solve(&a_neighborhood_geometry);
}

template <class CellType>
PlanarSeparator reconstructionWithLVIRA2D(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  LVIRA_2D<CellType> lvira_system;
  return lvira_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
}

template <class CellType>
PlanarSeparator reconstructionWithLVIRA3D(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  LVIRA_3D<CellType> lvira_system;
  return lvira_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
}

template <class CellType>
PlanarSeparator reconstructionWithMOF2D(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight, const double a_external_weight) {
  MOF_2D<CellType> mof_solver;
  return mof_solver.solve(
      CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>(&a_cell,
                                                                    &a_svm),
      a_internal_weight, a_external_weight);
}

template <class CellType>
PlanarSeparator reconstructionWithMOF3D(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight, const double a_external_weight) {
  MOF_3D<CellType> mof_solver;
  return mof_solver.solve(
      CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>(&a_cell,
                                                                    &a_svm),
      a_internal_weight, a_external_weight);
}

template <class MomentsContainerType, class CellType>
PlanarSeparator reconstructionWithAdvectedNormals(
    const MomentsContainerType& a_volume_moments_list,
    const R2PNeighborhood<CellType>& a_neighborhood,
    const double a_two_plane_threshold) {
  return AdvectedPlaneReconstruction::solve(
      a_volume_moments_list, a_neighborhood, a_two_plane_threshold);
}

//**********************************************************************
//     Function template definitions/inlined functions below this.
//     For debug versions.
//**********************************************************************

template <class CellType>
PlanarSeparator reconstructionWithR2P2DDebug(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  cleanReconstruction(
      a_neighborhood_geometry.getCenterCell(),
      (a_neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          a_neighborhood_geometry.getCenterCell().calculateVolume(),
      &a_initial_reconstruction);
  if (a_initial_reconstruction.getNumberOfPlanes() == 1) {
    R2PDebug<R2P_2D1P<CellType>> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  } else {
    R2PDebug<R2P_2D2P<CellType>> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  }
}

template <class CellType>
PlanarSeparator reconstructionWithR2P3DDebug(
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  cleanReconstruction(
      a_neighborhood_geometry.getCenterCell(),
      (a_neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          a_neighborhood_geometry.getCenterCell().calculateVolume(),
      &a_initial_reconstruction);
  if (a_initial_reconstruction.getNumberOfPlanes() == 1) {
    R2PDebug<R2P_3D1P<CellType>> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  } else {
    R2PDebug<R2P_3D2P<CellType>> r2p_system;
    return r2p_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
  }
}

PlanarSeparator reconstructionWithELVIRA2DDebug(
    const ELVIRANeighborhood& a_neighborhood_geometry) {
  ELVIRADebug<ELVIRA_2D> elvira_system;
  return elvira_system.solve(&a_neighborhood_geometry);
}

PlanarSeparator reconstructionWithELVIRA3DDebug(
    const ELVIRANeighborhood& a_neighborhood_geometry) {
  ELVIRADebug<ELVIRA_3D> elvira_system;
  return elvira_system.solve(&a_neighborhood_geometry);
}

template <class CellType>
PlanarSeparator reconstructionWithLVIRA2DDebug(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  LVIRADebug<LVIRA_2D<CellType>> lvira_system;
  return lvira_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
}

template <class CellType>
PlanarSeparator reconstructionWithLVIRA3DDebug(
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    PlanarSeparator a_initial_reconstruction) {
  LVIRADebug<LVIRA_3D<CellType>> lvira_system;
  return lvira_system.solve(a_neighborhood_geometry, a_initial_reconstruction);
}

template <class CellType>
PlanarSeparator reconstructionWithMOF2DDebug(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight, const double a_external_weight) {
  MOFDebug<MOF_2D<CellType>> mof_solver;
  return mof_solver.solve(
      CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>(&a_cell,
                                                                    &a_svm),
      a_internal_weight, a_external_weight);
}

template <class CellType>
PlanarSeparator reconstructionWithMOF3DDebug(
    const CellType& a_cell, const SeparatedMoments<VolumeMoments>& a_svm,
    const double a_internal_weight, const double a_external_weight) {
  MOFDebug<MOF_3D<CellType>> mof_solver;
  return mof_solver.solve(
      CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>(&a_cell,
                                                                    &a_svm),
      a_internal_weight, a_external_weight);
}

template <class MomentsContainerType, class CellType>
PlanarSeparator reconstructionWithAdvectedNormalsDebug(
    const MomentsContainerType& a_volume_moments_list,
    const R2PNeighborhood<CellType>& a_neighborhood,
    const double a_two_plane_threshold) {
  return AdvectedPlaneReconstructionDebug::solve(
      a_volume_moments_list, a_neighborhood, a_two_plane_threshold);
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_INTERFACE_TPP_
