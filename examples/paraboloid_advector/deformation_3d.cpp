// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/deformation_3d.h"

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

constexpr int NX = 256;
constexpr int NY = 256;
constexpr int NZ = 256;
constexpr int GC = 2;
constexpr IRL::Pt lower_domain(0.0, 0.0, 0.0);
constexpr IRL::Pt upper_domain(1.0, 1.0, 1.0);

BasicMesh Deformation3D::setMesh(void) {
  BasicMesh mesh(NX, NY, NZ, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Deformation3D::initialize(Data<double>* a_U, Data<double>* a_V,
                               Data<double>* a_W,
                               Data<IRL::Paraboloid>* a_interface) {
  Deformation3D::setVelocity(0.0, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt circle_center(1.0 * 0.35, 1.0 * 0.35, 1.0 * 0.35);
  const double circle_radius = 1.0 * 0.15;
  // const IRL::Pt circle_center(1.0 * 0.5, 1.0 * 0.5, 1.0 * 0.5);
  // const double circle_radius = 1.0 * 0.3;

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

void Deformation3D::setVelocity(const double a_time, Data<double>* a_U,
                                Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_U)(i, j, k) =
            2.0 * 1.0 * std::pow(sin(M_PI * mesh.xm(i) / 1.0), 2) *
            sin(2.0 * M_PI * mesh.ym(j) / 1.0) *
            sin(2.0 * M_PI * mesh.zm(k) / 1.0) * cos(M_PI * (a_time) / 3.0);
        (*a_V)(i, j, k) = -1.0 * std::pow(sin(M_PI * mesh.ym(j) / 1.0), 2) *
                          sin(2.0 * M_PI * mesh.xm(i) / 1.0) *
                          sin(2.0 * M_PI * mesh.zm(k) / 1.0) *
                          cos(M_PI * (a_time) / 3.0);
        (*a_W)(i, j, k) = -1.0 * std::pow(sin(M_PI * mesh.zm(k) / 1.0), 2) *
                          sin(2.0 * M_PI * mesh.xm(i) / 1.0) *
                          sin(2.0 * M_PI * mesh.ym(j) / 1.0) *
                          cos(M_PI * (a_time) / 3.0);
        // (*a_U)(i, j, k) = 1.0;
        // (*a_V)(i, j, k) = 1.0;
        // (*a_W)(i, j, k) = 1.0;
      }
    }
  }
}
