// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_SURFACE_OUTPUT_INTERSECTION_H_
#define IRL_GENERIC_CUTTING_PARABOLOID_SURFACE_OUTPUT_INTERSECTION_H_

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/ellipse.h"
#include "irl/paraboloid_reconstruction/parametrized_surface.h"

namespace IRL {
void addEllipseToSurfaceOutput(const AlignedParaboloid& a_aligned_paraboloid,
                               const Plane& a_face_plane,
                               const ParametrizedSurfaceOutput* a_surface);
}  // namespace IRL

#include "irl/generic_cutting/paraboloid_intersection/surface_output.tpp"

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_SURFACE_OUTPUT_INTERSECTION_H_
