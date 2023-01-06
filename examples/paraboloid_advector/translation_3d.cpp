// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/translation_3d.h"

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
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/moments/volume_moments.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/localized_separator_link.h"

#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/reconstruction_types.h"
#include "examples/paraboloid_advector/solver.h"
#include "examples/paraboloid_advector/vof_advection.h"

constexpr int NX = 10;
constexpr int NY = 10;
constexpr int NZ = 10;
constexpr int GC = 2;
constexpr IRL::Pt lower_domain(0.0, 0.0, 0.0);
constexpr IRL::Pt upper_domain(1.0, 1.0, 1.0);

BasicMesh Translation3D::setMesh(void) {
  BasicMesh mesh(NX, NY, NZ, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Translation3D::initialize(Data<double>* a_U, Data<double>* a_V,
                               Data<double>* a_W,
                               Data<IRL::Paraboloid>* a_interface) {
  Translation3D::setVelocity(0.0, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt circle_center(0.5 + 0.0 * mesh.dx(), 0.5 + 0.0 * mesh.dx(),
                              0.5 + 0.0 * mesh.dx());
  const double circle_radius = 0.25;

  // Loop over cells in domain. Skip if cell is not mixed phase.
  int sub_div = 10;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
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
          auto paraboloid =
              details::fromSphere(circle_center, circle_radius, circle_normal);

          double sub_dx =
              (mesh.x(i + 1) - mesh.x(i)) / static_cast<double>(sub_div);
          double sub_dy =
              (mesh.y(j + 1) - mesh.y(j)) / static_cast<double>(sub_div);
          double sub_dz =
              (mesh.z(k + 1) - mesh.z(k)) / static_cast<double>(sub_div);
          auto volume = IRL::Volume::fromScalarConstant(0.0);
          for (int kk = 0; kk < sub_div; ++kk) {
            for (int jj = 0; jj < sub_div; ++jj) {
              for (int ii = 0; ii < sub_div; ++ii) {
                const IRL::Pt lower_cell_pt(
                    mesh.x(i) + static_cast<double>(ii) * sub_dx,
                    mesh.y(j) + static_cast<double>(jj) * sub_dy,
                    mesh.z(k) + static_cast<double>(kk) * sub_dz);
                const IRL::Pt upper_cell_pt(
                    mesh.x(i) + static_cast<double>(ii + 1) * sub_dx,
                    mesh.y(j) + static_cast<double>(jj + 1) * sub_dy,
                    mesh.z(k) + static_cast<double>(kk + 1) * sub_dz);
                const auto sub_cell = IRL::RectangularCuboid::fromBoundingPts(
                    lower_cell_pt, upper_cell_pt);
                IRL::Normal normal =
                    sub_cell.calculateCentroid() - circle_center;
                normal.normalize();
                int largest_dir = 0;
                if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
                  largest_dir = 1;
                if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
                  largest_dir = 2;
                IRL::ReferenceFrame frame;
                if (largest_dir == 0)
                  frame[0] = crossProduct(normal, IRL::Normal(0.0, 1.0, 0.0));
                else if (largest_dir == 0)
                  frame[0] = crossProduct(normal, IRL::Normal(0.0, 0.0, 1.0));
                else
                  frame[0] = crossProduct(normal, IRL::Normal(1.0, 0.0, 0.0));
                frame[0].normalize();
                frame[1] = IRL::crossProduct(normal, frame[0]);
                frame[2] = normal;
                IRL::Paraboloid sub_paraboloid(
                    IRL::Pt(circle_center + circle_radius * normal), frame,
                    0.5 / circle_radius, 0.5 / circle_radius);
                volume += IRL::getVolumeMoments<IRL::Volume>(sub_cell,
                                                             sub_paraboloid);
              }
            }
          }
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, volume / cell.calculateVolume(), 1.0e-13,
                              paraboloid);
          auto new_frame = paraboloid.getReferenceFrame();
          auto new_datum =
              IRL::Pt(paraboloid.getDatum() +
                      solver_distance.getDistance() * new_frame[2]);
          paraboloid.setDatum(new_datum);
          (*a_interface)(i, j, k) = paraboloid;
        }
      }
    }
  }
  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void Translation3D::setVelocity(const double a_time, Data<double>* a_U,
                                Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_U)(i, j, k) = 1.0;
        (*a_V)(i, j, k) = 1.0 / 1.5;
        (*a_W)(i, j, k) = 1.0 / 3.0;
      }
    }
  }
}
