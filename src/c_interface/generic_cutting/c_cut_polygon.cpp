// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/generic_cutting/c_cut_polygon.h"

#include <cassert>

extern "C" {

void c_getPoly_RectCub_Poly(const c_RectCub* a_rectangular_cuboid,
                            const c_PlanarSep* a_separator,
                            const int* a_plane_index, c_Poly* a_polygon) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_polygon->obj_ptr) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
      *a_rectangular_cuboid->obj_ptr, *a_separator->obj_ptr,
      (*a_separator
            ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
}

void c_getPoly_Tet_Poly(const c_Tet* a_tet,
                            const c_PlanarSep* a_separator,
                            const int* a_plane_index, c_Poly* a_polygon) {
  assert(a_tet != nullptr);
  assert(a_tet->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_polygon->obj_ptr) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
      *a_tet->obj_ptr, *a_separator->obj_ptr,
      (*a_separator
            ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
}

void c_getPoly_Hex_Poly(const c_Hex* a_hexahedron,
                            const c_PlanarSep* a_separator,
                            const int* a_plane_index, c_Poly* a_polygon) {
  assert(a_hexahedron != nullptr);
  assert(a_hexahedron->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_polygon->obj_ptr) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
      *a_hexahedron->obj_ptr, *a_separator->obj_ptr,
      (*a_separator
            ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
}

void c_getPoly_RectCub_DivPoly(const c_RectCub* a_rectangular_cuboid,
                               const c_PlanarSep* a_separator,
                               const int* a_plane_index,
                               c_DivPoly* a_divided_polygon) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_divided_polygon != nullptr);
  assert(a_divided_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_divided_polygon->obj_ptr) =
      IRL::getPlanePolygonFromReconstruction<IRL::DividedPolygon>(
          *a_rectangular_cuboid->obj_ptr, *a_separator->obj_ptr,
          (*a_separator
                ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
  a_divided_polygon->obj_ptr->resetCentroid();
}

double c_getSA_RectCub(const c_RectCub* a_rectangular_cuboid,
                       const c_PlanarSep* a_separator) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  return IRL::getReconstructionSurfaceArea(*a_rectangular_cuboid->obj_ptr,
                                           *a_separator->obj_ptr);
}

}  // end extern C
