// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/interface_reconstruction_methods/c_volume_fraction_matching.h"
#include "src/data_structures/raw_pointer_wrapper.h"

extern "C" {

void c_matchVolumeFraction_RectCub_PlanarSep_Default(
    const c_RectCub* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr);
}

void c_matchVolumeFraction_RectCub_PlanarSep(
    const c_RectCub* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  assert(*a_volume_fraction_tolerance > 0.0);
  assert(*a_volume_fraction_tolerance < 1.0);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr,
                                        *a_volume_fraction_tolerance);
}

void c_matchVolumeFraction_Hex_PlanarSep_Default(
    const c_Hex* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr);
}

void c_matchVolumeFraction_Hex_PlanarSep(
    const c_Hex* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  assert(*a_volume_fraction_tolerance > 0.0);
  assert(*a_volume_fraction_tolerance < 1.0);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr,
                                        *a_volume_fraction_tolerance);
}

void c_matchVolumeFraction_Tet_PlanarSep_Default(
    const c_Tet* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr);
}

void c_matchVolumeFraction_Tet_PlanarSep(
    const c_Tet* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  assert(*a_volume_fraction_tolerance > 0.0);
  assert(*a_volume_fraction_tolerance < 1.0);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr,
                                        *a_volume_fraction_tolerance);
}

void c_matchVolumeFraction_TriPrism_PlanarSep_Default(
    const c_TriPrism* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr);
}

void c_matchVolumeFraction_TriPrism_PlanarSep(
    const c_TriPrism* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  assert(*a_volume_fraction_tolerance > 0.0);
  assert(*a_volume_fraction_tolerance < 1.0);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr,
                                        *a_volume_fraction_tolerance);
}

void c_matchVolumeFraction_Pyrmd_PlanarSep_Default(
    const c_Pyrmd* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr);
}

void c_matchVolumeFraction_Pyrmd_PlanarSep(
    const c_Pyrmd* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance) {
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->obj_ptr != nullptr);
  assert(*a_volume_fraction_tolerance > 0.0);
  assert(*a_volume_fraction_tolerance < 1.0);
  IRL::setDistanceToMatchVolumeFraction(*a_cell->obj_ptr, *a_volume_fraction,
                                        a_reconstruction->obj_ptr,
                                        *a_volume_fraction_tolerance);
}

void c_matchGroupVolumeFraction_Tet_PlanarSepPathGroup_Default(
    const c_Tet* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr,  wrapped_vf, a_reconstruction->obj_ptr);
}

void c_matchGroupVolumeFraction_Tet_PlanarSepPathGroup(
    const c_Tet* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	assert(*a_volume_fraction_tolerance > 0.0);
	assert(*a_volume_fraction_tolerance < 1.0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr, wrapped_vf, a_reconstruction->obj_ptr, *a_volume_fraction_tolerance);
}

void c_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_Default(
    const c_Pyrmd* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr,  wrapped_vf, a_reconstruction->obj_ptr);
}

void c_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup(
    const c_Pyrmd* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	assert(*a_volume_fraction_tolerance > 0.0);
	assert(*a_volume_fraction_tolerance < 1.0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr, wrapped_vf, a_reconstruction->obj_ptr, *a_volume_fraction_tolerance);
}

void c_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_Default(
    const c_TriPrism* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr,  wrapped_vf, a_reconstruction->obj_ptr);
}

void c_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup(
    const c_TriPrism* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	assert(*a_volume_fraction_tolerance > 0.0);
	assert(*a_volume_fraction_tolerance < 1.0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr, wrapped_vf, a_reconstruction->obj_ptr, *a_volume_fraction_tolerance);
}

void c_matchGroupVolumeFraction_Hex_PlanarSepPathGroup_Default(
    const c_Hex* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr,  wrapped_vf, a_reconstruction->obj_ptr);
}

void c_matchGroupVolumeFraction_Hex_PlanarSepPathGroup(
    const c_Hex* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance){
	assert(a_cell != nullptr);
	assert(a_cell->obj_ptr != nullptr);
	assert(a_reconstruction != nullptr);
	assert(a_reconstruction->obj_ptr != nullptr);
	assert(*a_size >= 0);
	assert(*a_volume_fraction_tolerance > 0.0);
	assert(*a_volume_fraction_tolerance < 1.0);
	auto wrapped_vf = IRL::wrapRawPointer(a_volume_fraction, static_cast<IRL::UnsignedIndex_t>(*a_size));
	IRL::setGroupDistanceToMatchVolumeFraction(*a_cell->obj_ptr, wrapped_vf, a_reconstruction->obj_ptr, *a_volume_fraction_tolerance);
}




}  // end extern C
