// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/interface_reconstruction_methods/c_reconstruction_interface.h"

#include <cassert>

extern "C" {

void c_reconstructELVIRA2D(const c_ELVIRANeigh* a_elvira_neighborhood,
                           c_PlanarSep* a_separator) {
  assert(a_elvira_neighborhood != nullptr);
  assert(a_elvira_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr =
      reconstructionWithELVIRA2D(*a_elvira_neighborhood->obj_ptr);
}

void c_reconstructELVIRA3D(const c_ELVIRANeigh* a_elvira_neighborhood,
                           c_PlanarSep* a_separator) {
  assert(a_elvira_neighborhood != nullptr);
  assert(a_elvira_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr =
      reconstructionWithELVIRA3D(*a_elvira_neighborhood->obj_ptr);
}

void c_reconstructMOF2D_RectCub(const c_RectCub* a_cell,
                                const c_SepVM* a_separated_volume_moments,
                                c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF2D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr);
}

void c_reconstructMOF3D_RectCub(const c_RectCub* a_cell,
                                const c_SepVM* a_separated_volume_moments,
                                c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF3D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr);
}

void c_reconstructMOF2D_GW_RectCub(const c_RectCub* a_cell,
                                   const c_SepVM* a_separated_volume_moments,
                                   const double* a_internal_weight,
                                   const double* a_external_weight,
                                   c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF2D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr,
      *a_internal_weight, *a_external_weight);
}

void c_reconstructMOF3D_GW_RectCub(const c_RectCub* a_cell,
                                   const c_SepVM* a_separated_volume_moments,
                                   const double* a_internal_weight,
                                   const double* a_external_weight,
                                   c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF3D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr,
      *a_internal_weight, *a_external_weight);
}

void c_reconstructMOF3D_Hex(const c_Hex* a_cell,
                                const c_SepVM* a_separated_volume_moments,
                                c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF3D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr);
}

void c_reconstructMOF3D_GW_Hex(const c_Hex* a_cell,
                                   const c_SepVM* a_separated_volume_moments,
                                   const double* a_internal_weight,
                                   const double* a_external_weight,
                                   c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF3D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr,
      *a_internal_weight, *a_external_weight);
}

void c_reconstructMOF2D_Tri(const c_Tri* a_cell,
                            const c_SepVM* a_separated_volume_moments,
                            c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF2D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr);
}

void c_reconstructMOF2D_GW_Tri(const c_Tri* a_cell,
                               const c_SepVM* a_separated_volume_moments,
                               const double* a_internal_weight,
                               const double* a_external_weight,
                               c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF2D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr,
      *a_internal_weight, *a_external_weight);
}

void c_reconstructMOF3D_Tet(const c_Tet* a_cell,
                            const c_SepVM* a_separated_volume_moments,
                            c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF3D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr);
}

void c_reconstructMOF3D_GW_Tet(const c_Tet* a_cell,
                               const c_SepVM* a_separated_volume_moments,
                               const double* a_internal_weight,
                               const double* a_external_weight,
                               c_PlanarSep* a_separator) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  *a_separator->obj_ptr = reconstructionWithMOF3D(
      *a_cell->obj_ptr, *a_separated_volume_moments->obj_ptr,
      *a_internal_weight, *a_external_weight);
}

void c_reconstructAdvectedNormals_RectCub(
    const c_ListVM_VMAN* a_volume_moments_list,
    const c_R2PNeigh_RectCub* a_neighborhood,
    const double* a_two_plane_threshold, c_PlanarSep* a_separator) {
  assert(a_volume_moments_list != nullptr);
  assert(a_volume_moments_list->obj_ptr != nullptr);
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithAdvectedNormals(
      *a_volume_moments_list->obj_ptr, *a_neighborhood->obj_ptr,
      *a_two_plane_threshold);
}

void c_reconstructAdvectedNormalsDbg_RectCub(
    const c_ListVM_VMAN* a_volume_moments_list,
    const c_R2PNeigh_RectCub* a_neighborhood,
    const double* a_two_plane_threshold, c_PlanarSep* a_separator) {
  assert(a_volume_moments_list != nullptr);
  assert(a_volume_moments_list->obj_ptr != nullptr);
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithAdvectedNormalsDebug(
      *a_volume_moments_list->obj_ptr, *a_neighborhood->obj_ptr,
      *a_two_plane_threshold);
}

void c_reconstructR2P2D_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithR2P2D(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructR2P3D_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithR2P3D(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructR2P2DDbg_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                   c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithR2P2DDebug(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructR2P3DDbg_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                   c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithR2P3DDebug(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructLVIRA2D_RectCub(const c_LVIRANeigh_RectCub* a_neighborhood,
                                  c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithLVIRA2D(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructLVIRA3D_RectCub(const c_LVIRANeigh_RectCub* a_neighborhood,
                                  c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithLVIRA3D(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructLVIRA2D_Hex(const c_LVIRANeigh_Hex* a_neighborhood,
                                  c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithLVIRA2D(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructLVIRA3D_Hex(const c_LVIRANeigh_Hex* a_neighborhood,
                                  c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithLVIRA3D(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}

void c_reconstructLVIRA3D_Tet(const c_LVIRANeigh_Tet* a_neighborhood,
                                  c_PlanarSep* a_separator) {
  assert(a_neighborhood != nullptr);
  assert(a_neighborhood->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  (*a_separator->obj_ptr) = IRL::reconstructionWithLVIRA3D(
      *a_neighborhood->obj_ptr, *a_separator->obj_ptr);
}


}  // end extern C
