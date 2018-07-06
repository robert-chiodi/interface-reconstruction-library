// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GENERIC_CUTTING_C_GENERIC_CUTTING_H_
#define SRC_C_INTERFACE_GENERIC_CUTTING_C_GENERIC_CUTTING_H_

#include "src/c_interface/geometry/polygons/c_polygon.h"
#include "src/c_interface/geometry/polygons/c_tri.h"
#include "src/c_interface/geometry/polyhedrons/c_capped_dodecahedron.h"
#include "src/c_interface/geometry/polyhedrons/c_capped_dodecahedron_doubles3.h"
#include "src/c_interface/geometry/polyhedrons/c_dodecahedron.h"
#include "src/c_interface/geometry/polyhedrons/c_hexahedron.h"
#include "src/c_interface/geometry/polyhedrons/c_octahedron.h"
#include "src/c_interface/geometry/polyhedrons/c_polyhedron24.h"
#include "src/c_interface/geometry/polyhedrons/c_polyhedron24_doubles3.h"
#include "src/c_interface/geometry/polyhedrons/c_pyramid.h"
#include "src/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "src/c_interface/geometry/polyhedrons/c_tet.h"
#include "src/c_interface/geometry/polyhedrons/c_triangular_prism.h"
#include "src/c_interface/geometry/polyhedrons/capped_dodecahedron_variations/c_capped_dodecahedron_LLLL.h"
#include "src/c_interface/geometry/polyhedrons/capped_dodecahedron_variations/c_capped_dodecahedron_LLLT.h"
#include "src/c_interface/geometry/polyhedrons/capped_dodecahedron_variations/c_capped_dodecahedron_LLTT.h"
#include "src/c_interface/geometry/polyhedrons/capped_dodecahedron_variations/c_capped_dodecahedron_LTLT.h"
#include "src/c_interface/geometry/polyhedrons/capped_dodecahedron_variations/c_capped_dodecahedron_LTTT.h"
#include "src/c_interface/geometry/polyhedrons/capped_dodecahedron_variations/c_capped_dodecahedron_TTTT.h"
#include "src/c_interface/geometry/polyhedrons/capped_octahedron_variations/c_capped_octahedron_LLL.h"
#include "src/c_interface/geometry/polyhedrons/capped_octahedron_variations/c_capped_octahedron_LLT.h"
#include "src/c_interface/geometry/polyhedrons/capped_octahedron_variations/c_capped_octahedron_LTT.h"
#include "src/c_interface/geometry/polyhedrons/capped_octahedron_variations/c_capped_octahedron_TTT.h"
#include "src/c_interface/geometry/polyhedrons/symmetric_decompositions/c_symmetric_hexahedron.h"
#include "src/c_interface/geometry/polyhedrons/symmetric_decompositions/c_symmetric_pyramid.h"
#include "src/c_interface/geometry/polyhedrons/symmetric_decompositions/c_symmetric_tet.h"
#include "src/c_interface/geometry/polyhedrons/symmetric_decompositions/c_symmetric_triangular_prism.h"
#include "src/c_interface/moments/c_separated_volume.h"
#include "src/c_interface/moments/c_separated_volume_moments.h"
#include "src/c_interface/moments/c_separated_volume_moments_doubles3.h"
#include "src/c_interface/moments/c_tagged_accumulated_listed_volume_moments_and_normal.h"
#include "src/c_interface/moments/c_tagged_accumulated_separated_volume.h"
#include "src/c_interface/moments/c_tagged_accumulated_separated_volume_moments.h"
#include "src/c_interface/moments/c_tagged_accumulated_tagged_accumulated_volume.h"
#include "src/c_interface/moments/c_tagged_accumulated_volume.h"
#include "src/c_interface/moments/c_tagged_accumulated_volume_moments.h"
#include "src/c_interface/planar_reconstruction/c_localized_separator.h"
#include "src/c_interface/planar_reconstruction/c_localized_separator_group_link.h"
#include "src/c_interface/planar_reconstruction/c_localized_separator_link.h"
#include "src/c_interface/planar_reconstruction/c_localizer_link.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/generic_cutting/generic_cutting.h"

extern "C" {

/// \file c_generic_cutting.h
///
/// These C-style funcions are
/// mapped to functions available in src/generic_cutting.h.
///
/// This file deals with functions that compute volume moments
/// for polyhedra and subdivided polyhedra.  In principle, the
/// first argument to the function is a pointer to a known
/// polytope class available in IRL, such as a
/// Polygon, Tet, or a Dodecahedron. The second argument
/// is a pointer to a PlanarSeparator, PlanarLocalizer,
/// LocalizedSeparator, or LocalizedSeparatorLink that
/// will subdivide or otherwise restrict the integration
/// area when calculating the volumetric moments. The third argument
/// is a pointer to an object of the  type of VolumeMoments that will be
/// returned.
///
/// Individual documentation for each
/// function is given alongside the function.

/// \brief Function to set the method used for cutting when a
/// c_getNormMoments function is called.
/// - 0 : RecursiveSimplexCutting
/// - 1 : HalfEdgeCutting
/// - 2 : SimplexCutting
void c_getMoments_setMethod(const int* a_cutting_method);

void c_getNormMoments_Dod_LocSepLink_SepVM(
    const c_Dod* a_Dod, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return);

void c_getNormMoments_Dod_LocSepLink_Vol(
    const c_Dod* a_Dod, const c_LocSepLink* a_localized_separator_link,
    double* a_moments_to_return);

void c_getNormMoments_TriPrism_LocSepLink_Vol(
    const c_TriPrism* a_triangular_prism,
    const c_LocSepLink* a_localized_separator_link,
    double* a_moments_to_return);

void c_getNormMoments_Octa_LocSepLink_Vol(
    const c_Octa* a_octahedron, const c_LocSepLink* a_localized_separator_link,
    double* a_moments_to_return);

void c_getNormMoments_CapDod_LocSepLink_SepVM(
    const c_CapDod* a_Cap_Dod, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return);

void c_getNormMoments_CapDod_d3_LocSepLink_SepVM_d3(
    const c_CapDod_d3* a_Cap_Dod,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM_d3* a_moments_to_return);

void c_getNormMoments_Poly24_LocSepLink_SepVM(
    const c_Poly24* a_polyhedron_24,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return);

void c_getNormMoments_Poly24_d3_LocSepLink_SepVM_d3(
    const c_Poly24_d3* a_polyhedron_24,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM_d3* a_moments_to_return);

void c_getMoments_CapDod_LocSepLink_SepVM(
    const c_CapDod* a_Cap_Dod, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return);

void c_getMoments_Dod_LocSepLink_SepVM(
    const c_Dod* a_Dod, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return);

void c_getMoments_Poly24_LocSepLink_SepVM(
    const c_Poly24* a_polyhedron_24,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return);

void c_getNormMoments_Tet_LocSepLink_SepVM(
    const c_Tet* a_tet, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return);

void c_getNormMoments_RectCub_PlanarSep_Vol(
    const c_RectCub* a_rectangular_cuboid,
    const c_PlanarSep* a_planar_separator, double* a_moments_to_return);

void c_getNormMoments_Tet_PlanarSep_Vol(const c_Tet* a_tet,
                                        const c_PlanarSep* a_planar_separator,
                                        double* a_moments_to_return);

void c_getNormMoments_TriPrism_PlanarSep_Vol(
    const c_TriPrism* a_tri_prism, const c_PlanarSep* a_planar_separator,
    double* a_moments_to_return);

void c_getNormMoments_Pyrmd_PlanarSep_Vol(const c_Pyrmd* a_pyramid,
                                          const c_PlanarSep* a_planar_separator,
                                          double* a_moments_to_return);

void c_getNormMoments_Hex_PlanarSep_Vol(const c_Hex* a_hexahedron,
                                        const c_PlanarSep* a_planar_separator,
                                        double* a_moments_to_return);

void c_getNormMoments_Hex_PlanarSep_SepVM(const c_Hex* a_hexahedron,
                                          const c_PlanarSep* a_planar_separator,
                                          c_SepVM* a_moments_to_return);

void c_getNormMoments_Dod_PlanarSep_SepVM(const c_Dod* a_Dod,
                                          const c_PlanarSep* a_planar_separator,
                                          c_SepVM* a_moments_to_return);

void c_getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM(
    const c_CapDod* a_Cap_Dod, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getNormMoments_Dod_LocSepLink_TagAccVM_SepVM(
    const c_Dod* a_Dod, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getMoments_Dod_LocSepLink_TagAccVM_SepVM(
    const c_Dod* a_Dod, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getNormMoments_Dod_LocSepLink_TagAccVM_SepVol(
    const c_Dod* a_Dod, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_Tet_LocSepLink_TagAccVM_SepVol(
    const c_Tet* a_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getMoments_Octa_LocSepLink_TagAccVM_SepVM(
    const c_Octa* a_octahedron, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getNormMoments_Octa_LocSepLink_TagAccVM_SepVol(
    const c_Octa* a_octahedron, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_RectCub_PlanarSep_SepVM(
    const c_RectCub* a_rectangular_cuboid,
    const c_PlanarSep* a_planar_separator, c_SepVM* a_moments_to_return);

void c_getNormMoments_Tri_LocLink_TagAccVM_VM(
    const c_Tri* a_tri, const c_LocLink* a_localizer_link,
    c_TagAccVM_VM* a_moments_to_return);

void c_getNormMoments_Tri_PlanarLoc_Vol(const c_Tri* a_tri,
                                        const c_PlanarLoc* a_planar_localizer,
                                        double* a_moments_to_return);

void c_getNormMoments_Poly_PlanarLoc_Vol(const c_Poly* a_poly,
                                         const c_PlanarLoc* a_planar_localizer,
                                         double* a_moments_to_return);

void c_getMoments_Tri_LocLink_TagAccListVM_VMAN(
    const c_Tri* a_tri, const c_LocLink* a_localizer_link,
    c_TagAccListVM_VMAN* a_moments_to_return);

void c_getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_LLL* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_LLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_LTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_TTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymTet* a_sym_tet,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymHex* a_sym_hex,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol(
    const c_SymTet* a_sym_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol(
    const c_SymHex* a_sym_hex, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return);

void c_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM(
    const c_SymTet* a_sym_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM(
    const c_SymHex* a_sym_hex, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getMoments_SymTet_LocSepLink_TagAccVM_SepVM(
    const c_SymTet* a_sym_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getMoments_SymHex_LocSepLink_TagAccVM_SepVM(
    const c_SymHex* a_sym_hex, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return);

void c_getNormMoments_CapDod_LLLL_LocSepLink_SepVol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLT_LocSepLink_SepVol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LTLT_LocSepLink_SepVol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLTT_LocSepLink_SepVol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LTTT_LocSepLink_SepVol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_TTTT_LocSepLink_SepVol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLL_LocSepLink_SepVol(
    const c_CapOcta_LLL* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLT_LocSepLink_SepVol(
    const c_CapOcta_LLT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LTT_LocSepLink_SepVol(
    const c_CapOcta_LTT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_TTT_LocSepLink_SepVol(
    const c_CapOcta_TTT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_Tet_LocSepLink_SepVol(
    const c_Tet* a_tet, const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_LLL* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_LLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_LTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_TTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_LLL* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_LLT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_LTT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_TTT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol(
    const c_Tet* a_tet, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol(
    const c_Pyrmd* a_pyramid, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol(
    const c_TriPrism* a_tri_prism, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol(
    const c_Hex* a_hex, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return);

void c_getNormMoments_SymTet_LocSep_SepVol(
    const c_SymTet* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_SymPyrmd_LocSep_SepVol(
    const c_SymPyrmd* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_SymTriPrism_LocSep_SepVol(
    const c_SymTriPrism* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_SymHex_LocSep_SepVol(
    const c_SymHex* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLL_LocSep_SepVol(
    const c_CapOcta_LLL* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LLT_LocSep_SepVol(
    const c_CapOcta_LLT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_LTT_LocSep_SepVol(
    const c_CapOcta_LTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapOcta_TTT_LocSep_SepVol(
    const c_CapOcta_TTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLL_LocSep_SepVol(
    const c_CapDod_LLLL* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLLT_LocSep_SepVol(
    const c_CapDod_LLLT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LTLT_LocSep_SepVol(
    const c_CapDod_LTLT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LLTT_LocSep_SepVol(
    const c_CapDod_LLTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_LTTT_LocSep_SepVol(
    const c_CapDod_LTTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

void c_getNormMoments_CapDod_TTTT_LocSep_SepVol(
    const c_CapDod_TTTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return);

}  // end extern C

#endif  // SRC_C_INTERFACE_GENERIC_CUTTING_C_GENERIC_CUTTING_H_
