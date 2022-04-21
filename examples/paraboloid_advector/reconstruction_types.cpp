// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/reconstruction_types.h"

#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/parameters/constants.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/vof_advection.h"

void getReconstruction(
    const std::string& a_reconstruction_method,
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedParaboloidLink>& a_localized_paraboloid_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::Paraboloid>* a_interface) {
  if (a_reconstruction_method == "KnownCircle") {
    KnownCircle::getReconstruction(a_liquid_volume_fraction, a_dt, a_U, a_V,
                                   a_W, a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: KnownCircle. \n";
    std::exit(-1);
  }
}

void KnownCircle::getReconstruction(
    const Data<double>& a_liquid_volume_fraction, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::Paraboloid>* a_interface) {
  const BasicMesh& mesh = a_U.getMesh();
  static constexpr double circle_radius = 0.15;
  static IRL::Pt circle_center(0.0, 0.25, mesh.zm(mesh.kmin()));
  // First order Forward Euler advection of circle center
  circle_center = back_project_vertex(circle_center, a_dt, a_U, a_V, a_W);

  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (a_liquid_volume_fraction(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          continue;
        } else {
          // Add known reconstruction here for now
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
          auto circle_normal =
              IRL::Normal::fromPtNormalized(mid_pt - circle_center);
          circle_normal[2] = 0.0;
          (*a_interface)(i, j, k) =
              details::fromCircle(circle_center, circle_radius, circle_normal);
        }
      }
    }
  }
  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void correctInterfacePlaneBorders(Data<IRL::Paraboloid>* a_interface) {
  const BasicMesh& mesh = (*a_interface).getMesh();
  // Fix distances in reconstruction for periodic boundary

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[0] -= mesh.lx();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[0] += mesh.lx();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[1] -= mesh.ly();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[1] += mesh.ly();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[2] -= mesh.lz();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        IRL::Pt datum = (*a_interface)(i, j, k).getDatum();
        datum[2] += mesh.lz();
        (*a_interface)(i, j, k).setDatum(datum);
      }
    }
  }
}
