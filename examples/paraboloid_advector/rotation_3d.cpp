// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/rotation_3d.h"

#include <float.h>
#include <chrono>
#include <cmath>
#include <iostream>

#include "irl/distributions/k_means.h"
#include "irl/distributions/partition_by_normal_vector.h"
#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/moments/volume_moments.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/localized_separator_link.h"

#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/reconstruction_types.h"
#include "examples/paraboloid_advector/solver.h"
#include "examples/paraboloid_advector/vof_advection.h"

constexpr int NX = 30;
constexpr int NY = 30;
constexpr int NZ = NX / 4 + 4;
constexpr int GC = 2;
constexpr IRL::Pt lower_domain(-0.5, -0.5,
                               -0.5 * static_cast<double>(NZ) /
                                   static_cast<double>(NX));
constexpr IRL::Pt upper_domain(0.5, 0.5,
                               0.5 * static_cast<double>(NZ) /
                                   static_cast<double>(NX));

BasicMesh Rotation3D::setMesh(void) {
  BasicMesh mesh(NX, NY, NZ, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Rotation3D::initialize(Data<double>* a_U, Data<double>* a_V,
                            Data<double>* a_W,
                            Data<IRL::Paraboloid>* a_interface) {
  Rotation3D::setVelocity(0.0, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt circle_center(0.0, 0.25, 0.0);
  const double circle_radius = 0.125;

  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const IRL::Pt mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
        IRL::Pt disp = mid_pt - circle_center;
        const auto mag = magnitude(disp);
        if (mag < circle_radius - 2.0 * mesh.dx()) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
        } else if (mag > circle_radius + 2.0 * mesh.dx()) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else {
          auto circle_normal = IRL::Normal::fromPt(disp);
          circle_normal.normalize();
          (*a_interface)(i, j, k) =
              details::fromSphere(circle_center, circle_radius, circle_normal);
        }
      }
    }
  }
  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void Rotation3D::setVelocity(const double a_time, Data<double>* a_U,
                             Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  auto rotation_axis = IRL::Normal(1.0, 1.0, 1.0);
  const double angular_velocity = 2.0 * M_PI;
  rotation_axis.normalize();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt position(mesh.xm(i) - 0.5, mesh.ym(j) - 0.5, mesh.zm(k) - 0.5);
        position -= IRL::dotProduct(rotation_axis, position) * rotation_axis;
        IRL::Normal velocity = IRL::crossProduct(rotation_axis, position);
        velocity.normalize();
        velocity *= angular_velocity * IRL::magnitude(position);
        // (*a_U)(i, j, k) = velocity[0];
        // (*a_V)(i, j, k) = velocity[1];
        // (*a_W)(i, j, k) = velocity[2];
        (*a_U)(i, j, k) = -2.0 * M_PI * mesh.ym(j);
        (*a_V)(i, j, k) = 2.0 * M_PI * mesh.xm(i);
        (*a_W)(i, j, k) = 0.0;
      }
    }
  }
}
