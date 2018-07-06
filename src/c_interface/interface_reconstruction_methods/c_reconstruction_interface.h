// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_RECONSTRUCTION_INTERFACE_H_
#define SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_RECONSTRUCTION_INTERFACE_H_

#include "src/c_interface/geometry/polygons/c_tri.h"
#include "src/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "src/c_interface/geometry/polyhedrons/c_hexahedron.h"
#include "src/c_interface/geometry/polyhedrons/c_tet.h"
#include "src/c_interface/interface_reconstruction_methods/c_elvira_neighborhood.h"
#include "src/c_interface/interface_reconstruction_methods/c_lvira_neighborhood_hexahedron.h"
#include "src/c_interface/interface_reconstruction_methods/c_lvira_neighborhood_tet.h"
#include "src/c_interface/interface_reconstruction_methods/c_lvira_neighborhood_rectangular_cuboid.h"
#include "src/c_interface/interface_reconstruction_methods/c_r2p_neighborhood_rectangular_cuboid.h"
#include "src/c_interface/moments/c_listedvm_vman.h"
#include "src/c_interface/moments/c_separated_volume_moments.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/geometry/polygons/tri.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/interface_reconstruction_methods/reconstruction_interface.h"
#include "src/planar_reconstruction/planar_separator.h"

extern "C" {
/// \file c_localizers.h
///
/// These C-style funcions are
/// mapped to functions available in src/reconstruction_interface.h.
///
/// This file includes functions to place PlanarSeparator objects
/// in geometries. These methods differ in what they require. For
/// the individual needs of each reconstruction method,
/// it is best to constult its specific documentation.

void c_reconstructELVIRA2D(const c_ELVIRANeigh* a_elvira_neighborhood,
                           c_PlanarSep* a_separator);

void c_reconstructELVIRA3D(const c_ELVIRANeigh* a_elvira_neighborhood,
                           c_PlanarSep* a_separator);

void c_reconstructMOF2D_RectCub(const c_RectCub* a_cell,
                                const c_SepVM* a_separated_volume_moments,
                                c_PlanarSep* a_separator);

void c_reconstructMOF3D_RectCub(const c_RectCub* a_cell,
                                const c_SepVM* a_separated_volume_moments,
                                c_PlanarSep* a_separator);

void c_reconstructMOF2D_GW_RectCub(const c_RectCub* a_cell,
                                   const c_SepVM* a_separated_volume_moments,
                                   const double* a_internal_weight,
                                   const double* a_external_weight,
                                   c_PlanarSep* a_separator);

void c_reconstructMOF3D_GW_RectCub(const c_RectCub* a_cell,
                                   const c_SepVM* a_separated_volume_moments,
                                   const double* a_internal_weight,
                                   const double* a_external_weight,
                                   c_PlanarSep* a_separator);

void c_reconstructMOF3D_Hex(const c_Hex* a_cell,
                            const c_SepVM* a_separated_volume_moments,
                            c_PlanarSep* a_separator);

void c_reconstructMOF3D_GW_Hex(const c_Hex* a_cell,
                               const c_SepVM* a_separated_volume_moments,
                               const double* a_internal_weight,
                               const double* a_external_weight,
                               c_PlanarSep* a_separator);

void c_reconstructMOF2D_Tri(const c_Tri* a_cell,
                            const c_SepVM* a_separated_volume_moments,
                            c_PlanarSep* a_separator);

void c_reconstructMOF2D_GW_Tri(const c_Tri* a_cell,
                               const c_SepVM* a_separated_volume_moments,
                               const double* a_internal_weight,
                               const double* a_external_weight,
                               c_PlanarSep* a_separator);

void c_reconstructMOF3D_Tet(const c_Tet* a_cell,
                            const c_SepVM* a_separated_volume_moments,
                            c_PlanarSep* a_separator);

void c_reconstructMOF3D_GW_Tet(const c_Tet* a_cell,
                               const c_SepVM* a_separated_volume_moments,
                               const double* a_internal_weight,
                               const double* a_external_weight,
                               c_PlanarSep* a_separator);

void c_reconstructAdvectedNormals_RectCub(
    const c_ListVM_VMAN* a_volume_moments_list,
    const c_R2PNeigh_RectCub* a_neighborhood,
    const double* a_two_plane_threshold, c_PlanarSep* a_separator);

void c_reconstructAdvectedNormalsDbg_RectCub(
    const c_ListVM_VMAN* a_volume_moments_list,
    const c_R2PNeigh_RectCub* a_neighborhood,
    const double* a_two_plane_threshold, c_PlanarSep* a_separator);

void c_reconstructR2P2D_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                c_PlanarSep* a_separator);

void c_reconstructR2P3D_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                c_PlanarSep* a_separator);

void c_reconstructR2P2DDbg_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                   c_PlanarSep* a_separator);

void c_reconstructR2P3DDbg_RectCub(const c_R2PNeigh_RectCub* a_neighborhood,
                                   c_PlanarSep* a_separator);

void c_reconstructLVIRA2D_RectCub(const c_LVIRANeigh_RectCub* a_neighborhood,
                                  c_PlanarSep* a_separator);

void c_reconstructLVIRA3D_RectCub(const c_LVIRANeigh_RectCub* a_neighborhood,
                                  c_PlanarSep* a_separator);

void c_reconstructLVIRA2D_Hex(const c_LVIRANeigh_Hex* a_neighborhood,
                                  c_PlanarSep* a_separator);

void c_reconstructLVIRA3D_Hex(const c_LVIRANeigh_Hex* a_neighborhood,
                                  c_PlanarSep* a_separator);

void c_reconstructLVIRA3D_Tet(const c_LVIRANeigh_Tet* a_neighborhood,
                                  c_PlanarSep* a_separator);

}  // end extern C

#endif  // SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_RECONSTRUCTION_INTERFACE_H_
