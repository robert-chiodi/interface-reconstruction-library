// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/advector/deformation_2d.h"

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

#include "examples/advector/data.h"
#include "examples/advector/reconstruction_types.h"
#include "examples/advector/solver.h"
#include "examples/advector/vof_advection.h"

constexpr int GC = 2;
constexpr IRL::Pt lower_domain(0.0, 0.0, 0.0);
constexpr IRL::Pt upper_domain(1.0, 1.0, 1.0);

BasicMesh Deformation2D::setMesh(const IRL::UnsignedIndex_t a_n) {
  BasicMesh mesh(a_n, a_n, 1, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  const double dx =
      (my_upper_domain[0] - my_lower_domain[0]) / static_cast<double>(a_n);
  my_lower_domain[2] = 0.0 - dx;
  my_upper_domain[2] = 0.0 + dx;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Deformation2D::initialize(Data<double>* a_U, Data<double>* a_V,
                               Data<double>* a_W,
                               Data<IRL::PlanarSeparator>* a_separators) {
  Deformation2D::setVelocity(0.0, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt circle_center(0.5, 0.75, mesh.zm(mesh.kmin()));
  const double circle_radius = 0.15;
  constexpr int subdivisions = 1;
  IRL::PlanarSeparator temp_separator;
  IRL::ListedVolumeMoments<IRL::VolumeMomentsAndNormal> listed_volume_moments;
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        IRL::RectangularCuboid cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        double dx = upper_cell_pt.x() - lower_cell_pt.x();
        double dy = upper_cell_pt.y() - lower_cell_pt.y();
        double sc_dx = dx / static_cast<double>(subdivisions);
        double sc_dy = dy / static_cast<double>(subdivisions);
        IRL::Pt sc_lower;
        IRL::Pt sc_upper;
        if (IRL::magnitude(cell.calculateCentroid() - circle_center) -
                circle_radius >
            3.0 * dx) {
          (*a_separators)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), -10000000.0));
        } else if (IRL::magnitude(cell.calculateCentroid() - circle_center) -
                       circle_radius <
                   -3.0 * dx) {
          (*a_separators)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), 10000000.0));
        } else {
          listed_volume_moments.clear();
          // Create separators and localizers for sub-divided cell
          for (int ii = 0; ii < subdivisions; ++ii) {
            for (int jj = 0; jj < subdivisions; ++jj) {
              sc_lower[0] = lower_cell_pt[0] + static_cast<double>(ii) * sc_dx;
              sc_lower[1] = lower_cell_pt[1] + static_cast<double>(jj) * sc_dy;
              sc_lower[2] = lower_cell_pt[2];
              sc_upper[0] = sc_lower[0] + sc_dx;
              sc_upper[1] = sc_lower[1] + sc_dy;
              sc_upper[2] = upper_cell_pt[2];
              IRL::RectangularCuboid sub_cell =
                  IRL::RectangularCuboid::fromBoundingPts(sc_lower, sc_upper);
              IRL::Normal sub_cell_normal = IRL::Normal::fromPtNormalized(
                  (sub_cell.calculateCentroid() - circle_center));
              temp_separator = IRL::PlanarSeparator::fromOnePlane(IRL::Plane(
                  sub_cell_normal,
                  sub_cell_normal * IRL::Pt(circle_center +
                                            IRL::Normal::toPt(sub_cell_normal *
                                                              circle_radius))));
              IRL::Polygon interface_poly =
                  IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                      sub_cell, temp_separator, temp_separator[0]);
              auto moments = IRL::getVolumeMoments<IRL::VolumeMomentsAndNormal>(
                  interface_poly);
              listed_volume_moments += moments;
            }
          }
          IRL::VolumeMomentsAndNormal mean_moments;

          for (const auto& element : listed_volume_moments) {
            mean_moments += element;
          }
          if (mean_moments.volumeMoments().volume() == 0.0) {
            (*a_separators)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
                IRL::Plane(IRL::Normal(0.0, 0.0, 0.0),
                           std::copysign(1.0, circle_radius -
                                                  IRL::magnitude(
                                                      cell.calculateCentroid() -
                                                      circle_center))));
          } else {
            mean_moments.normalizeByVolume();
            mean_moments.normal().normalize();
            (*a_separators)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
                IRL::Plane(mean_moments.normal(),
                           mean_moments.normal() *
                               mean_moments.volumeMoments().centroid()));
          }
        }
      }
    }
  }
  a_separators->updateBorder();
}

void Deformation2D::setVelocity(const double a_time, Data<double>* a_U,
                                Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_U)(i, j, k) = -2.0 * std::pow(sin(M_PI * mesh.xm(i)), 2) *
                          sin(M_PI * mesh.ym(j)) * cos(M_PI * mesh.ym(j)) *
                          cos(M_PI * (a_time) / 8.0);
        (*a_V)(i, j, k) = 2.0 * std::pow(sin(M_PI * mesh.ym(j)), 2) *
                          sin(M_PI * mesh.xm(i)) * cos(M_PI * mesh.xm(i)) *
                          cos(M_PI * (a_time) / 8.0);
        (*a_W)(i, j, k) = 0.0;
      }
    }
  }
}
