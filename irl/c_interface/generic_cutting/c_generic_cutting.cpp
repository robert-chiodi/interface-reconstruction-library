// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/generic_cutting/c_generic_cutting.h"

#include <cassert>

namespace IRL {

enum class c_RuntimeCuttingMethod {
  RecursiveSimplexCutting = 0,
  HalfEdgeCutting = 1,
  SimplexCutting = 2,
};

#ifdef C_STATIC_CUTTING
template <class ReturnType, class EncompassingType, class ReconstructionType>
static ReturnType c_RuntimegetNormMoments(
    const EncompassingType& a_encompassing_polytope,
    const ReconstructionType& a_reconstruction,
    const c_RuntimeCuttingMethod& a_cutting_method) {
  return IRL::getNormMoments<ReturnType, DefaultCuttingMethod>(
      a_encompassing_polytope, a_reconstruction);
}

template <class ReturnType, class EncompassingType, class ReconstructionType>
static ReturnType c_RuntimegetMoments(
    const EncompassingType& a_encompassing_polytope,
    const ReconstructionType& a_reconstruction,
    const c_RuntimeCuttingMethod& a_cutting_method) {
  return IRL::getVolumeMoments<ReturnType, DefaultCuttingMethod>(
      a_encompassing_polytope, a_reconstruction);
}

#else   // C_STATIC_CUTTING not defined
template <class ReturnType, class EncompassingType, class ReconstructionType>
static ReturnType c_RuntimegetNormMoments(
    const EncompassingType& a_encompassing_polytope,
    const ReconstructionType& a_reconstruction,
    const c_RuntimeCuttingMethod& a_cutting_method) {
  switch (a_cutting_method) {
    case c_RuntimeCuttingMethod::RecursiveSimplexCutting:
      return IRL::getNormalizedVolumeMoments<ReturnType,
                                             RecursiveSimplexCutting>(
          a_encompassing_polytope, a_reconstruction);

    case c_RuntimeCuttingMethod::HalfEdgeCutting:
      return IRL::getNormalizedVolumeMoments<ReturnType, HalfEdgeCutting>(
          a_encompassing_polytope, a_reconstruction);

    case c_RuntimeCuttingMethod::SimplexCutting:
      return IRL::getNormalizedVolumeMoments<ReturnType, SimplexCutting>(
          a_encompassing_polytope, a_reconstruction);
    default:
      std::cout << "During call to cutting: Unkown cutting method required for "
                   "getNormalizedVolumeMoments in "
                   "the C/Fortran interface."
                << std::endl;
      std::exit(-1);
  }
}

template <class ReturnType, class EncompassingType, class ReconstructionType>
static ReturnType c_RuntimegetMoments(
    const EncompassingType& a_encompassing_polytope,
    const ReconstructionType& a_reconstruction,
    const c_RuntimeCuttingMethod& a_cutting_method) {
  switch (a_cutting_method) {
    case c_RuntimeCuttingMethod::RecursiveSimplexCutting:
      return IRL::getVolumeMoments<ReturnType, RecursiveSimplexCutting>(
          a_encompassing_polytope, a_reconstruction);

    case c_RuntimeCuttingMethod::HalfEdgeCutting:
      return IRL::getVolumeMoments<ReturnType, HalfEdgeCutting>(
          a_encompassing_polytope, a_reconstruction);

    case c_RuntimeCuttingMethod::SimplexCutting:
      return IRL::getVolumeMoments<ReturnType, SimplexCutting>(
          a_encompassing_polytope, a_reconstruction);
    default:
      std::cout << "During call to cutting: Unkown cutting method required for "
                   "getNormalizedVolumeMoments in "
                   "the C/Fortran interface."
                << std::endl;
      std::exit(-1);
  }
}
#endif  // C_STATIC_CUTTING

template <class CuttingType>
constexpr c_RuntimeCuttingMethod c_getCompiledDefaultCuttingMethod(void);

template <>
constexpr c_RuntimeCuttingMethod
c_getCompiledDefaultCuttingMethod<RecursiveSimplexCutting>(void) {
  return c_RuntimeCuttingMethod::RecursiveSimplexCutting;
}

template <>
constexpr c_RuntimeCuttingMethod
c_getCompiledDefaultCuttingMethod<HalfEdgeCutting>(void) {
  return c_RuntimeCuttingMethod::HalfEdgeCutting;
}

template <>
constexpr c_RuntimeCuttingMethod
c_getCompiledDefaultCuttingMethod<SimplexCutting>(void) {
  return c_RuntimeCuttingMethod::SimplexCutting;
}

// Initialize C_CUTTING_METHOD with the default cutting method in
// src/generic_cutting/default_cutting_method.h
static IRL::c_RuntimeCuttingMethod C_CUTTING_METHOD =
    c_getCompiledDefaultCuttingMethod<DefaultCuttingMethod>();

}  // namespace IRL

extern "C" {

#ifdef C_STATIC_CUTTING
void c_getMoments_setMethod(const int* a_cutting_method) {
  std::cout << "Attempted to set a cutting method during runtime while the "
               "C/Fortran interface to the IRL library was compiled to use the "
               "DefaultCuttingMethod defined in "
               "irl/generic_cutting/default_cutting_method.h"
            << std::endl;
  std::exit(-1);
}
#else   // C_STATIC_CUTTING not defined
void c_getMoments_setMethod(const int* a_cutting_method) {
  switch (*a_cutting_method) {
    case 0:
      IRL::C_CUTTING_METHOD =
          IRL::c_RuntimeCuttingMethod::RecursiveSimplexCutting;
      break;
    case 1:
      IRL::C_CUTTING_METHOD = IRL::c_RuntimeCuttingMethod::HalfEdgeCutting;
      break;
    case 2:
      IRL::C_CUTTING_METHOD = IRL::c_RuntimeCuttingMethod::SimplexCutting;
      break;
    default:
      std::cout
          << "Unkown cutting method required for getNormalizedVolumeMoments in "
             "the C/Fortran interface."
          << std::endl;
      std::exit(-1);
  }
}
#endif  // C_STATIC_CUTTING

void c_getNormMoments_Dod_LocSepLink_SepVM(
    const c_Dod* a_dodecahedron, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return) {
  assert(a_dodecahedron != nullptr);
  assert(a_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Dod_LocSepLink_Vol(
    const c_Dod* a_dodecahedron, const c_LocSepLink* a_localized_separator_link,
    double* a_moments_to_return) {
  assert(a_dodecahedron != nullptr);
  assert(a_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD));
}

void c_getNormMoments_Octa_LocSepLink_Vol(
    const c_Octa* a_octahedron, const c_LocSepLink* a_localized_separator_link,
    double* a_moments_to_return) {
  assert(a_octahedron != nullptr);
  assert(a_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD));
}

void c_getNormMoments_TriPrism_LocSepLink_Vol(
    const c_TriPrism* a_triangular_prism,
    const c_LocSepLink* a_localized_separator_link,
    double* a_moments_to_return) {
  assert(a_triangular_prism != nullptr);
  assert(a_triangular_prism->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_triangular_prism->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD));
}

void c_getNormMoments_CapDod_LocSepLink_SepVM(
    const c_CapDod* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_d3_LocSepLink_SepVM_d3(
    const c_CapDod_d3* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM_d3* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetNormMoments<
      IRL::SeparatedMoments<IRL::VolumeMomentsAndDoubles<3>>>(
      *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Poly24_LocSepLink_SepVM(
    const c_Poly24* a_polyhedron_24,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return) {
  assert(a_polyhedron_24 != nullptr);
  assert(a_polyhedron_24->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_polyhedron_24->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Poly24_d3_LocSepLink_SepVM_d3(
    const c_Poly24_d3* a_polyhedron_24,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM_d3* a_moments_to_return) {
  assert(a_polyhedron_24 != nullptr);
  assert(a_polyhedron_24->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetNormMoments<
      IRL::SeparatedMoments<IRL::VolumeMomentsAndDoubles<3>>>(
      *a_polyhedron_24->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getMoments_CapDod_LocSepLink_SepVM(
    const c_CapDod* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_Dod_LocSepLink_SepVM(
    const c_Dod* a_dodecahedron, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return) {
  assert(a_dodecahedron != nullptr);
  assert(a_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_Poly24_LocSepLink_SepVM(
    const c_Poly24* a_polyhedron_24,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return) {
  assert(a_polyhedron_24 != nullptr);
  assert(a_polyhedron_24->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_polyhedron_24->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Tet_LocSepLink_SepVM(
    const c_Tet* a_tet, const c_LocSepLink* a_localized_separator_link,
    c_SepVM* a_moments_to_return) {
  assert(a_tet != nullptr);
  assert(a_tet->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_tet->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_RectCub_PlanarSep_Vol(
    const c_RectCub* a_rectangular_cuboid,
    const c_PlanarSep* a_planar_separator, double* a_moments_to_return) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_rectangular_cuboid->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD));
}

void c_getNormMoments_Tet_PlanarSep_Vol(const c_Tet* a_tet,
                                        const c_PlanarSep* a_planar_separator,
                                        double* a_moments_to_return) {
  assert(a_tet != nullptr);
  assert(a_tet->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_tet->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD));
}

void c_getNormMoments_TriPrism_PlanarSep_Vol(
    const c_TriPrism* a_tri_prism, const c_PlanarSep* a_planar_separator,
    double* a_moments_to_return) {
  assert(a_tri_prism != nullptr);
  assert(a_tri_prism->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_tri_prism->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD));
}
void c_getNormMoments_Pyrmd_PlanarSep_Vol(const c_Pyrmd* a_pyramid,
                                          const c_PlanarSep* a_planar_separator,
                                          double* a_moments_to_return) {
  assert(a_pyramid != nullptr);
  assert(a_pyramid->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_pyramid->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD));
}

void c_getNormMoments_Hex_PlanarSep_Vol(const c_Hex* a_hexahedron,
                                        const c_PlanarSep* a_planar_separator,
                                        double* a_moments_to_return) {
  assert(a_hexahedron != nullptr);
  assert(a_hexahedron->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return =
      static_cast<double>(IRL::c_RuntimegetNormMoments<IRL::Volume>(
          *a_hexahedron->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD));
}

void c_getNormMoments_Hex_PlanarSep_SepVM(const c_Hex* a_hexahedron,
                                          const c_PlanarSep* a_planar_separator,
                                          c_SepVM* a_moments_to_return) {
  assert(a_hexahedron != nullptr);
  assert(a_hexahedron->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_hexahedron->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Dod_PlanarSep_SepVM(const c_Dod* a_dodecahedron,
                                          const c_PlanarSep* a_planar_separator,
                                          c_SepVM* a_moments_to_return) {
  assert(a_dodecahedron != nullptr);
  assert(a_dodecahedron->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_dodecahedron->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LocSepLink_TagAccVM_SepVM(
    const c_CapDod* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_CapDod_LocSepLink_TagAccVM_SepVM(
    const c_CapDod* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Dod_LocSepLink_TagAccVM_SepVM(
    const c_Dod* a_dodecahedron, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_dodecahedron != nullptr);
  assert(a_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_Dod_LocSepLink_TagAccVM_SepVM(
    const c_Dod* a_dodecahedron, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_dodecahedron != nullptr);
  assert(a_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Dod_LocSepLink_TagAccVM_SepVol(
    const c_Dod* a_dodecahedron, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_dodecahedron != nullptr);
  assert(a_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Tet_LocSepLink_TagAccVM_SepVol(
    const c_Tet* a_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_tet != nullptr);
  assert(a_tet->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_tet->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getMoments_Octa_LocSepLink_TagAccVM_SepVM(
    const c_Octa* a_octahedron, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_octahedron != nullptr);
  assert(a_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Octa_LocSepLink_TagAccVM_SepVol(
    const c_Octa* a_octahedron, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_octahedron != nullptr);
  assert(a_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_RectCub_PlanarSep_SepVM(
    const c_RectCub* a_rectangular_cuboid,
    const c_PlanarSep* a_planar_separator, c_SepVM* a_moments_to_return) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
          *a_rectangular_cuboid->obj_ptr, *a_planar_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Tri_LocLink_TagAccVM_VM(
    const c_Tri* a_tri, const c_LocLink* a_localizer_link,
    c_TagAccVM_VM* a_moments_to_return) {
  assert(a_tri != nullptr);
  assert(a_tri->obj_ptr != nullptr);
  assert(a_localizer_link != nullptr);
  assert(a_localizer_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetNormMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::VolumeMoments>>(
      *a_tri->obj_ptr, *a_localizer_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Tri_PlanarLoc_Vol(const c_Tri* a_tri,
                                        const c_PlanarLoc* a_planar_localizer,
                                        double* a_moments_to_return) {
  assert(a_tri != nullptr);
  assert(a_tri->obj_ptr != nullptr);
  assert(a_planar_localizer != nullptr);
  assert(a_planar_localizer->obj_ptr != nullptr);
  (*a_moments_to_return) = IRL::c_RuntimegetMoments<IRL::Volume>(
      *a_tri->obj_ptr, *a_planar_localizer->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Poly_PlanarSep_Vol(const c_Poly* a_poly,
                                         const c_PlanarSep* a_planar_separator,
                                         double* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  (*a_moments_to_return) = IRL::c_RuntimegetMoments<IRL::Volume>(
      *a_poly->obj_ptr, *a_planar_separator->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Poly_PlanarLoc_Vol(const c_Poly* a_poly,
                                         const c_PlanarLoc* a_planar_localizer,
                                         double* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_planar_localizer != nullptr);
  assert(a_planar_localizer->obj_ptr != nullptr);
  (*a_moments_to_return) = IRL::c_RuntimegetMoments<IRL::Volume>(
      *a_poly->obj_ptr, *a_planar_localizer->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Tri_PlanarLoc_VM(const c_Tri* a_tri,
                                       const c_PlanarLoc* a_planar_localizer,
                                       c_VM* a_moments_to_return) {
  assert(a_tri != nullptr);
  assert(a_tri->obj_ptr != nullptr);
  assert(a_planar_localizer != nullptr);
  assert(a_planar_localizer->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<IRL::VolumeMoments>(
      *a_tri->obj_ptr, *a_planar_localizer->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Poly_PlanarLoc_VM(const c_Poly* a_poly,
                                        const c_PlanarLoc* a_planar_localizer,
                                        c_VM* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_planar_localizer != nullptr);
  assert(a_planar_localizer->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<IRL::VolumeMoments>(
      *a_poly->obj_ptr, *a_planar_localizer->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getMoments_Tri_LocLink_TagAccListVM_VMAN(
    const c_Tri* a_tri, const c_LocLink* a_localizer_link,
    c_TagAccListVM_VMAN* a_moments_to_return) {
  assert(a_tri != nullptr);
  assert(a_tri->obj_ptr != nullptr);
  assert(a_localizer_link != nullptr);
  assert(a_localizer_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  (*a_moments_to_return->obj_ptr) = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedListedVolumeMoments<IRL::VolumeMomentsAndNormal>>(
      *a_tri->obj_ptr, *a_localizer_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLL_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTLT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLTT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTTT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_TTTT_LocSepLink_TagAccVM_SepVol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLL_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_LLL* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLT_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_LLT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LTT_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_LTT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_TTT_LocSepLink_TagAccVM_SepVol(
    const c_CapOcta_TTT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTet_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymTet* a_sym_tet,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_sym_tet != nullptr);
  assert(a_sym_tet->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_sym_tet->obj_ptr, *a_localized_separator_group_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymPyrmd_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_sym_pyramid != nullptr);
  assert(a_sym_pyramid->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_sym_pyramid->obj_ptr, *a_localized_separator_group_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTriPrism_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_sym_tri_prism != nullptr);
  assert(a_sym_tri_prism->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_sym_tri_prism->obj_ptr, *a_localized_separator_group_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymHex_LocSepGroupLink_TagAccVM2_Vol(
    const c_SymHex* a_sym_hex,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_sym_hex != nullptr);
  assert(a_sym_hex->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_sym_hex->obj_ptr, *a_localized_separator_group_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVol(
    const c_SymTet* a_sym_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_sym_tet != nullptr);
  assert(a_sym_tet->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_sym_tet->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVol(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_sym_pyramid != nullptr);
  assert(a_sym_pyramid->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_sym_pyramid->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVol(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_sym_tri_prism != nullptr);
  assert(a_sym_tri_prism->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_sym_tri_prism->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVol(
    const c_SymHex* a_sym_hex, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVol* a_moments_to_return) {
  assert(a_sym_hex != nullptr);
  assert(a_sym_hex->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::SeparatedMoments<IRL::Volume>>>(
      *a_sym_hex->obj_ptr, *a_localized_separator_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTet_LocSepLink_TagAccVM_SepVM(
    const c_SymTet* a_sym_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_tet != nullptr);
  assert(a_sym_tet->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_tet->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_pyramid != nullptr);
  assert(a_sym_pyramid->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_pyramid->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_tri_prism != nullptr);
  assert(a_sym_tri_prism->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_tri_prism->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymHex_LocSepLink_TagAccVM_SepVM(
    const c_SymHex* a_sym_hex, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_hex != nullptr);
  assert(a_sym_hex->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_hex->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_SymTet_LocSepLink_TagAccVM_SepVM(
    const c_SymTet* a_sym_tet, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_tet != nullptr);
  assert(a_sym_tet->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_tet->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_SymPyrmd_LocSepLink_TagAccVM_SepVM(
    const c_SymPyrmd* a_sym_pyramid,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_pyramid != nullptr);
  assert(a_sym_pyramid->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_pyramid->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_SymTriPrism_LocSepLink_TagAccVM_SepVM(
    const c_SymTriPrism* a_sym_tri_prism,
    const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_tri_prism != nullptr);
  assert(a_sym_tri_prism->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_tri_prism->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getMoments_SymHex_LocSepLink_TagAccVM_SepVM(
    const c_SymHex* a_sym_hex, const c_LocSepLink* a_localized_separator_link,
    c_TagAccVM_SepVM* a_moments_to_return) {
  assert(a_sym_hex != nullptr);
  assert(a_sym_hex->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          *a_sym_hex->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLL_LocSepLink_SepVol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLT_LocSepLink_SepVol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTLT_LocSepLink_SepVol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLTT_LocSepLink_SepVol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTTT_LocSepLink_SepVol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_TTTT_LocSepLink_SepVol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_dodecahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLL_LocSepLink_SepVol(
    const c_CapOcta_LLL* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLT_LocSepLink_SepVol(
    const c_CapOcta_LLT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LTT_LocSepLink_SepVol(
    const c_CapOcta_LTT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_TTT_LocSepLink_SepVol(
    const c_CapOcta_TTT* a_capped_octahedron,
    const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_capped_octahedron->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Tet_LocSepLink_SepVol(
    const c_Tet* a_tet, const c_LocSepLink* a_localized_separator_link,
    c_SepVol* a_moments_to_return) {
  assert(a_tet != nullptr);
  assert(a_tet->obj_ptr != nullptr);
  assert(a_localized_separator_link != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_tet->obj_ptr, *a_localized_separator_link->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_dodecahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_dodecahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_dodecahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_dodecahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_dodecahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_dodecahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_LLL* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_octahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_LLT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_octahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_LTT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_octahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM2_Vol(
    const c_CapOcta_TTT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM2_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
          *a_capped_octahedron->obj_ptr,
          *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLL_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LLLL* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr,
      *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LLLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr,
      *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTLT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LTLT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr,
      *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LLTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr,
      *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_LTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr,
      *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_TTTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapDod_TTTT* a_capped_dodecahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_dodecahedron != nullptr);
  assert(a_capped_dodecahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_dodecahedron->obj_ptr,
      *a_localized_separator_group_link->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLL_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_LLL* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_group_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_LLT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_group_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_LTT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_group_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_TTT_LocSepGroupLink_TagAccVM_Vol(
    const c_CapOcta_TTT* a_capped_octahedron,
    const c_LocSepGroupLink* a_localized_separator_group_link,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_capped_octahedron != nullptr);
  assert(a_capped_octahedron->obj_ptr != nullptr);
  assert(a_localized_separator_group_link != nullptr);
  assert(a_localized_separator_group_link->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::AccumulateWrapper<IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>>(
      *a_capped_octahedron->obj_ptr, *a_localized_separator_group_link->obj_ptr,
      IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Tet_PlanarSepPathGroup_TagAccVM_Vol(
    const c_Tet* a_tet, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_tet != nullptr);
  assert(a_tet->obj_ptr != nullptr);
  assert(a_path != nullptr);
  assert(a_path->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
      *a_tet->obj_ptr, *a_path->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Pyrmd_PlanarSepPathGroup_TagAccVM_Vol(
    const c_Pyrmd* a_pyramid, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_pyramid != nullptr);
  assert(a_pyramid->obj_ptr != nullptr);
  assert(a_path != nullptr);
  assert(a_path->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
      *a_pyramid->obj_ptr, *a_path->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_TriPrism_PlanarSepPathGroup_TagAccVM_Vol(
    const c_TriPrism* a_tri_prism, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_tri_prism != nullptr);
  assert(a_tri_prism->obj_ptr != nullptr);
  assert(a_path != nullptr);
  assert(a_path->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
      *a_tri_prism->obj_ptr, *a_path->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_Hex_PlanarSepPathGroup_TagAccVM_Vol(
    const c_Hex* a_hex, const c_PlanarSepPathGroup* a_path,
    c_TagAccVM_Vol* a_moments_to_return) {
  assert(a_hex != nullptr);
  assert(a_hex->obj_ptr != nullptr);
  assert(a_path != nullptr);
  assert(a_path->obj_ptr != nullptr);
  assert(a_moments_to_return != nullptr);
  assert(a_moments_to_return->obj_ptr != nullptr);
  *a_moments_to_return->obj_ptr = IRL::c_RuntimegetMoments<
      IRL::TaggedAccumulatedVolumeMoments<IRL::Volume>>(
      *a_hex->obj_ptr, *a_path->obj_ptr, IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTet_LocSep_SepVol(
    const c_SymTet* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymPyrmd_LocSep_SepVol(
    const c_SymPyrmd* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymTriPrism_LocSep_SepVol(
    const c_SymTriPrism* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_SymHex_LocSep_SepVol(
    const c_SymHex* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLL_LocSep_SepVol(
    const c_CapOcta_LLL* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LLT_LocSep_SepVol(
    const c_CapOcta_LLT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_LTT_LocSep_SepVol(
    const c_CapOcta_LTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapOcta_TTT_LocSep_SepVol(
    const c_CapOcta_TTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLL_LocSep_SepVol(
    const c_CapDod_LLLL* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLLT_LocSep_SepVol(
    const c_CapDod_LLLT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTLT_LocSep_SepVol(
    const c_CapDod_LTLT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LLTT_LocSep_SepVol(
    const c_CapDod_LLTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_LTTT_LocSep_SepVol(
    const c_CapDod_LTTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

void c_getNormMoments_CapDod_TTTT_LocSep_SepVol(
    const c_CapDod_TTTT* a_poly, const c_LocSep* a_localized_separator,
    c_SepVol* a_moments_to_return) {
  assert(a_poly != nullptr);
  assert(a_poly->obj_ptr != nullptr);
  assert(a_localized_separator != nullptr);
  assert(a_moments_to_return != nullptr);
  *a_moments_to_return->obj_ptr =
      IRL::c_RuntimegetNormMoments<IRL::SeparatedMoments<IRL::Volume>>(
          *a_poly->obj_ptr, *a_localized_separator->obj_ptr,
          IRL::C_CUTTING_METHOD);
}

}  // end extern C
