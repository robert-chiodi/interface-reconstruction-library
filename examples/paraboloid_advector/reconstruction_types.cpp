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
#include "irl/interface_reconstruction_methods/constrained_optimization_behavior.h"
#include "irl/interface_reconstruction_methods/elvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/plvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/moments/volume_moments_with_gradient.h"
#include "irl/moments/volume_with_gradient.h"
#include "irl/optimization/constrained_levenberg_marquardt.h"
#include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
#include "irl/paraboloid_reconstruction/hessian_paraboloid.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/planar_separator.h"

#include <Eigen/Dense>
#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/vof_advection.h"

void getReconstruction(const std::string& a_reconstruction_method,
                       const Data<double>& a_liquid_volume_fraction,
                       const Data<IRL::Pt>& a_liquid_centroid,
                       const Data<IRL::Pt>& a_gas_centroid,
                       const Data<IRL::LocalizedParaboloidLink<double>>&
                           a_localized_paraboloid_link,
                       const double a_dt, const Data<double>& a_U,
                       const Data<double>& a_V, const Data<double>& a_W,
                       Data<IRL::Paraboloid>* a_interface) {
  if (a_reconstruction_method == "Jibben") {
    Jibben::getReconstruction(a_liquid_volume_fraction, a_dt, a_U, a_V, a_W,
                              a_interface);
  } else if (a_reconstruction_method == "CentroidFit") {
    Centroid::getReconstruction(a_liquid_volume_fraction, a_dt, a_U, a_V, a_W,
                                a_interface);
  } else if (a_reconstruction_method == "PLIC") {
    PLIC::getReconstruction(a_liquid_volume_fraction, a_dt, a_U, a_V, a_W,
                            a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Valid entries are: PLIC, CentroidFit, Jibben. \n";
    std::exit(-1);
  }
}

// Wendland radial basis function
// Wendland, H. (1995). Piecewise polynomial, positive definite and
// compactly supported radial functions of minimal degree. Advances in
// Computational Mathematics, 4(1), 389â€“396.
double wgauss(const double d, const double h) {
  if (d >= h) {
    return 0.0;
  } else {
    return (1.0 + 4.0 * d / h) * std::pow(1.0 - d / h, 4.0);
  }
}

void updateReconstructionELVIRA(
    const Data<double>& a_liquid_volume_fraction,
    Data<IRL::PlanarSeparator>* a_liquid_gas_interface) {
  IRL::ELVIRANeighborhood neighborhood;
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  neighborhood.resize(27);
  IRL::RectangularCuboid cells[27];
  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
            a_liquid_volume_fraction(i, j, k) >
                IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, a_liquid_volume_fraction(i, j, k) - 0.5);
          (*a_liquid_gas_interface)(i, j, k) =
              IRL::PlanarSeparator::fromOnePlane(
                  IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }
        // Build surrounding stencil information for ELVIRA.
        for (int kk = k - 1; kk < k + 2; ++kk) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            for (int ii = i - 1; ii < i + 2; ++ii) {
              // Reversed order, bad for cache locality but thats okay..
              cells[(kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1)] =
                  IRL::RectangularCuboid::fromBoundingPts(
                      IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                      IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              neighborhood.setMember(
                  &cells[(kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1)],
                  &a_liquid_volume_fraction(ii, jj, kk), ii - i, jj - j,
                  kk - k);
            }
          }
        }
        // Now perform actual ELVIRA and obtain interface PlanarSeparator
        (*a_liquid_gas_interface)(i, j, k) =
            reconstructionWithELVIRA3D(neighborhood);
      }
    }
  }
  // Update border with simple ghost-cell fill and correct distances for
  // assumed periodic boundary
  a_liquid_gas_interface->updateBorder();
  // correctInterfacePlaneBorders(a_liquid_gas_interface);
}

// Reconstruction with LVIRA - use input PlanarSeparator as initial guess
void updateReconstructionLVIRA(
    const Data<double>& a_liquid_volume_fraction, const int a_nneigh,
    Data<IRL::PlanarSeparator>* a_liquid_gas_interface) {
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();

  IRL::LVIRANeighborhood<IRL::RectangularCuboid> neighborhood;
  std::vector<IRL::RectangularCuboid> cells;
  // std::vector<double> weights; // maybe later

  const int grid_size =
      (a_nneigh * 2 + 1) * (a_nneigh * 2 + 1) * (a_nneigh * 2 + 1);
  neighborhood.resize(grid_size);
  cells.resize(grid_size);

  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
            a_liquid_volume_fraction(i, j, k) >
                IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, a_liquid_volume_fraction(i, j, k) - 0.5);
          (*a_liquid_gas_interface)(i, j, k) =
              IRL::PlanarSeparator::fromOnePlane(
                  IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }

        // Build surrounding stencil information for LVIRA.
        IRL::UnsignedIndex_t ndata = 0;
        for (int kk = k - a_nneigh; kk <= k + a_nneigh; ++kk) {
          for (int jj = j - a_nneigh; jj <= j + a_nneigh; ++jj) {
            for (int ii = i - a_nneigh; ii <= i + a_nneigh; ++ii) {
              // Trap center cell
              if (ii == i && jj == j && kk == k) {
                neighborhood.setCenterOfStencil(ndata);
              }
              cells[ndata] = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              neighborhood.setMember(ndata, &cells[ndata],
                                     &a_liquid_volume_fraction(ii, jj, kk));
              // Increment counter
              ++ndata;
            }
          }
        }
        auto found_planar_separator = (*a_liquid_gas_interface)(i, j, k);
        // Now perform actual LVIRA and obtain interface PlanarSeparator
        (*a_liquid_gas_interface)(i, j, k) =
            reconstructionWithLVIRA3D(neighborhood, found_planar_separator);
      }
    }
  }
}

void updatePolygon(const Data<double>& a_liquid_volume_fraction,
                   const Data<IRL::PlanarSeparator>& a_liquid_gas_interface,
                   Data<IRL::Polygon>* a_interface_polygon) {
  const BasicMesh& mesh = a_liquid_gas_interface.getMesh();
  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
            a_liquid_volume_fraction(i, j, k) >
                IRL::global_constants::VF_HIGH) {
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        (*a_interface_polygon)(i, j, k) =
            IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                cell, a_liquid_gas_interface(i, j, k),
                a_liquid_gas_interface(i, j, k)[0]);
      }
    }
  }
}

std::array<double, 6> fitParaboloidToPLICHeights(
    const Data<IRL::Polygon>& a_polygon, const Data<double>& a_volume_fraction,
    const IRL::Pt& a_reference_point, const IRL::ReferenceFrame& a_frame,
    const int a_i, const int a_j, const int a_k, const int a_nneigh,
    const double a_width) {
  const BasicMesh& mesh = a_polygon.getMesh();
  const double meshsize = 1.0;  // mesh.dx();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
  const int ic(a_i), jc(a_j), kc(a_k);
  const IRL::Pt pref = a_reference_point;
  const auto frame = a_frame;

  for (int k = kc - a_nneigh; k <= kc + a_nneigh; ++k) {
    for (int j = jc - a_nneigh; j <= jc + a_nneigh; ++j) {
      for (int i = ic - a_nneigh; i <= ic + a_nneigh; ++i) {
        const IRL::UnsignedIndex_t shape =
            a_polygon(i, j, k).getNumberOfVertices();
        if (shape == 0) {
          continue;
        }
        // Local polygon normal and centroid
        IRL::Pt ploc = a_polygon(i, j, k).calculateCentroid();
        IRL::Normal nloc = a_polygon(i, j, k).calculateNormal();
        // if (frame[2] * nloc <= 0.0) {
        //   continue;
        // }
        ploc -= pref;
        ploc /= meshsize;
        const IRL::Pt tmp_pt = ploc;
        const IRL::Normal tmp_n = nloc;
        for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
          ploc[d] = frame[d] * tmp_pt;
          nloc[d] = frame[d] * tmp_n;
        }
        // Plane coefficients
        Eigen::VectorXd reconstruction_plane_coeffs(3);
        reconstruction_plane_coeffs << -(ploc * nloc) / meshsize, nloc[0],
            nloc[1];
        reconstruction_plane_coeffs /= -nloc[2];
        // Integrals
        Eigen::VectorXd integrals = Eigen::VectorXd::Zero(6);
        double b_dot_sum = 0.0;
        for (IRL::UnsignedIndex_t v = 0; v < shape; ++v) {
          IRL::UnsignedIndex_t vn = (v + 1) % shape;
          IRL::Pt vert1 = a_polygon(i, j, k)[v];
          IRL::Pt vert2 = a_polygon(i, j, k)[vn];
          vert1 -= pref;
          vert2 -= pref;
          vert1 /= meshsize;
          vert2 /= meshsize;
          IRL::Pt tmp_pt1 = vert1;
          IRL::Pt tmp_pt2 = vert2;
          for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
            vert1[d] = frame[d] * tmp_pt1;
            vert2[d] = frame[d] * tmp_pt2;
          }

          const double xv = vert1[0];
          const double yv = vert1[1];
          const double xvn = vert2[0];
          const double yvn = vert2[1];

          Eigen::VectorXd integral_to_add(6);
          integral_to_add << (xv * yvn - xvn * yv) / 2.0,
              (xv + xvn) * (xv * yvn - xvn * yv) / 6.0,
              (yv + yvn) * (xv * yvn - xvn * yv) / 6.0,
              (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0,
              (yvn - yv) *
                  (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
                   2.0 * xv * xvn * yvn + xvn * xvn * yv +
                   3.0 * xvn * xvn * yvn) /
                  24.0,
              (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
          integrals += integral_to_add;
        }
        b_dot_sum += integrals.head(3).dot(reconstruction_plane_coeffs);

        // Get weighting
        const double gaussianweight =  // 1.0;
            a_width <= 0.0
                ? 1.0
                : wgauss(std::sqrt(static_cast<IRL::Vec3<double>>(ploc) *
                                   static_cast<IRL::Vec3<double>>(ploc)),
                         a_width);
        const double vfrac = a_volume_fraction(i, j, k);
        double vfrac_weight = 1.0;
        const double limit_vfrac = 0.1;
        if (vfrac < limit_vfrac) {
          vfrac_weight = 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
        } else if (vfrac > 1.0 - limit_vfrac) {
          vfrac_weight =
              0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
        }
        double ww = 1.0;
        ww *= gaussianweight;
        ww *= vfrac_weight;

        if (ww > 0.0) {
          A += ww * integrals * integrals.transpose();
          b += ww * integrals * b_dot_sum;
        }
      }
    }
  }
  Eigen::VectorXd sol = A.colPivHouseholderQr().solve(b);
  return std::array<double, 6>{
      {sol(0), sol(1), sol(2), sol(3), sol(4), sol(5)}};
}

std::array<double, 6> fitParaboloidToCentroids(
    const Data<IRL::Polygon>& a_polygon, const Data<double>& a_volume_fraction,
    const IRL::Pt& a_reference_point, const IRL::ReferenceFrame& a_frame,
    const int a_i, const int a_j, const int a_k, const int a_nneigh,
    const double a_width) {
  const BasicMesh& mesh = a_polygon.getMesh();
  const double meshsize = 1.0;  // mesh.dx();
  // const int ncells = std::pow(2 * a_nneigh + 1, 3);
  const int ncells =
      (2 * a_nneigh + 1) * (2 * a_nneigh + 1) * (2 * a_nneigh + 1);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(ncells, 6);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(ncells);
  const int ic(a_i), jc(a_j), kc(a_k);
  const IRL::Pt pref = a_reference_point;
  const auto frame = a_frame;

  IRL::UnsignedIndex_t ndata = 0;
  for (int k = kc - a_nneigh; k <= kc + a_nneigh; ++k) {
    for (int j = jc - a_nneigh; j <= jc + a_nneigh; ++j) {
      for (int i = ic - a_nneigh; i <= ic + a_nneigh; ++i) {
        if (a_polygon(i, j, k).getNumberOfVertices() == 0) {
          continue;
        }
        IRL::Pt ploc = a_polygon(i, j, k).calculateCentroid();
        const IRL::Normal nloc = a_polygon(i, j, k).calculateNormal();
        const double surf = a_polygon(i, j, k).calculateAbsoluteVolume() /
                            (meshsize * meshsize);
        const double normalproj = std::max(frame[2] * nloc, 0.0);
        // if (normalproj <= 0.0) continue;
        ploc -= pref;
        ploc /= meshsize;
        const IRL::Pt tmp_pt = ploc;
        for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
          ploc[d] = frame[d] * tmp_pt;
        }
        const double gaussianweight =
            // 1.0;
            a_width <= 0.0
                ? 1.0
                : wgauss(std::sqrt(static_cast<IRL::Vec3<double>>(ploc) *
                                   static_cast<IRL::Vec3<double>>(ploc)),
                         a_width);

        const double vfrac = a_volume_fraction(i, j, k);
        double vfrac_weight = 1.0;
        const double limit_vfrac = 0.1;
        if (vfrac < limit_vfrac) {
          vfrac_weight = 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
        } else if (vfrac > 1.0 - limit_vfrac) {
          vfrac_weight =
              0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
        }
        double ww = 1.0;
        ww *= normalproj;
        ww *= surf;
        ww *= gaussianweight;
        ww *= vfrac_weight;

        if (ww > 0.0) {
          // Store least squares matrix and RHS
          A(ndata, 0) = std::sqrt(ww);
          A(ndata, 1) = 0.0;
          A(ndata, 2) = 0.0;
          A(ndata, 3) = std::sqrt(ww) * ploc[0] * ploc[0];
          A(ndata, 4) = std::sqrt(ww) * ploc[0] * ploc[1];
          A(ndata, 5) = std::sqrt(ww) * ploc[1] * ploc[1];
          b(ndata) = std::sqrt(ww) * ploc[2];
          // Increment counter
          ++ndata;
        }
      }
    }
  }
  A.conservativeResize(ndata, Eigen::NoChange);
  b.conservativeResize(ndata, Eigen::NoChange);
  Eigen::VectorXd sol = A.colPivHouseholderQr().solve(b);
  return std::array<double, 6>{
      // {sol(0), sol(1), sol(2), sol(3), sol(4), sol(5)}};
      {sol(0), 0.0, 0.0, sol(3), sol(4), sol(5)}};
}

void Jibben::getReconstruction(const Data<double>& a_liquid_volume_fraction,
                               const double a_dt, const Data<double>& a_U,
                               const Data<double>& a_V, const Data<double>& a_W,
                               Data<IRL::Paraboloid>* a_interface) {
  const BasicMesh& mesh = a_U.getMesh();

  Data<IRL::PlanarSeparator> interface(&mesh);
  updateReconstructionELVIRA(a_liquid_volume_fraction, &interface);
  updateReconstructionLVIRA(a_liquid_volume_fraction, 1, &interface);
  Data<IRL::Polygon> polygon(&mesh);
  updatePolygon(a_liquid_volume_fraction, interface, &polygon);
  polygon.updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] -= mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] += mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] -= mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] += mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] -= mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] += mesh.lz();
        }
      }
    }
  }

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (a_liquid_volume_fraction(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          // continue;
        } else {
          const IRL::Normal norm_poly = polygon(i, j, k).calculateNormal();
          const double poly_area = polygon(i, j, k).calculateVolume();
          const IRL::Pt pref = polygon(i, j, k).calculateCentroid();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[1]))
            largest_dir = 1;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(norm_poly, fit_frame[0]);
          fit_frame[2] = norm_poly;
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Paraboloid paraboloid;

          double sum_vfrac = 0.0;
          for (int kk = -1; kk < 2; ++kk) {
            for (int jj = -1; jj < 2; ++jj) {
              for (int ii = -1; ii < 2; ++ii) {
                sum_vfrac += a_liquid_volume_fraction(i + ii, j + jj, k + kk);
              }
            }
          }

          auto sol_fit =
              fitParaboloidToPLICHeights(polygon, a_liquid_volume_fraction,
                                         pref, fit_frame, i, j, k, 1, 0.0);
          const double a = sol_fit[0], b = sol_fit[1], c = sol_fit[2],
                       d = sol_fit[3], e = sol_fit[4], f = sol_fit[5];
          const double theta = 0.5 * std::atan2(e, (IRL::safelyTiny(d - f)));
          const double cos_t = std::cos(theta);
          const double sin_t = std::sin(theta);
          const double A =
              -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
          const double B =
              -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
          // Translation to coordinate system R' where aligned paraboloid
          // valid Translation is R' = {x' = x + u, y' = y + v, z' = z + w}
          const double denominator = IRL::safelyTiny(4.0 * d * f - e * e);
          const double u = (2.0 * b * f - c * e) / denominator;
          const double v = -(b * e - 2.0 * d * c) / denominator;
          const double w =
              -(a + (-b * b * f + b * c * e - c * c * d) / denominator);

          IRL::UnitQuaternion rotation(theta, fit_frame[2]);
          IRL::Pt datum =
              pref - u * fit_frame[0] - v * fit_frame[1] - w * fit_frame[2];
          auto new_frame = rotation * fit_frame;
          const double max_curvature_dx = 1.0;
          double a_coeff = A;
          double b_coeff = B;
          if (std::sqrt(u * u + v * v + w * w) > 10.0 * mesh.dx() ||
              std::fabs(A) * mesh.dx() > max_curvature_dx ||
              std::fabs(B) * mesh.dx() > max_curvature_dx) {
            paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
          } else {
            if (fabs(a_coeff) < 1.0e-3) {
              a_coeff = std::copysign(1.0e-3, a_coeff);
            }
            if (fabs(b_coeff) < 1.0e-3) {
              b_coeff = std::copysign(1.0e-3, b_coeff);
            }
            paraboloid = IRL::Paraboloid(datum, new_frame, a_coeff, b_coeff);
          }

          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, a_liquid_volume_fraction(i, j, k), 1.0e-14,
                              paraboloid);

          if (solver_distance.getDistance() == -DBL_MAX) {
            paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
            IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
                new_solver_distance(cell, a_liquid_volume_fraction(i, j, k),
                                    1.0e-14, paraboloid);
            if (new_solver_distance.getDistance() == -DBL_MAX) {
              (*a_interface)(i, j, k) =
                  IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
            } else {
              auto new_datum =
                  IRL::Pt(paraboloid.getDatum() +
                          new_solver_distance.getDistance() * fit_frame[2]);
              paraboloid.setDatum(new_datum);
              (*a_interface)(i, j, k) = paraboloid;
            }
          } else {
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            (*a_interface)(i, j, k) = paraboloid;
          }
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void Centroid::getReconstruction(const Data<double>& a_liquid_volume_fraction,
                                 const double a_dt, const Data<double>& a_U,
                                 const Data<double>& a_V,
                                 const Data<double>& a_W,
                                 Data<IRL::Paraboloid>* a_interface) {
  const BasicMesh& mesh = a_U.getMesh();

  Data<IRL::PlanarSeparator> interface(&mesh);
  updateReconstructionELVIRA(a_liquid_volume_fraction, &interface);
  updateReconstructionLVIRA(a_liquid_volume_fraction, 1, &interface);
  Data<IRL::Polygon> polygon(&mesh);
  updatePolygon(a_liquid_volume_fraction, interface, &polygon);
  polygon.updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] -= mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] += mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] -= mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] += mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] -= mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] += mesh.lz();
        }
      }
    }
  }

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (a_liquid_volume_fraction(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          // continue;
        } else {
          const IRL::Normal norm_poly = polygon(i, j, k).calculateNormal();
          const double poly_area = polygon(i, j, k).calculateVolume();
          const IRL::Pt pref = polygon(i, j, k).calculateCentroid();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[1]))
            largest_dir = 1;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(norm_poly, fit_frame[0]);
          fit_frame[2] = norm_poly;
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Paraboloid paraboloid;

          double sum_vfrac = 0.0;
          for (int kk = -1; kk < 2; ++kk) {
            for (int jj = -1; jj < 2; ++jj) {
              for (int ii = -1; ii < 2; ++ii) {
                sum_vfrac += a_liquid_volume_fraction(i + ii, j + jj, k + kk);
              }
            }
          }

          if (std::fabs(sum_vfrac - 27.0 * 0.5) < 10.0 &&
              a_liquid_volume_fraction(i, j, k) > 1.0e-6 &&
              a_liquid_volume_fraction(i, j, k) < 1.0 - 1.0e-6) {
            auto sol_fit =
                fitParaboloidToCentroids(polygon, a_liquid_volume_fraction,
                                         pref, fit_frame, i, j, k, 1, 0.0);
            const double a = sol_fit[0], b = sol_fit[1], c = sol_fit[2],
                         d = sol_fit[3], e = sol_fit[4], f = sol_fit[5];
            const double theta = 0.5 * std::atan2(e, (IRL::safelyTiny(d - f)));
            const double cos_t = std::cos(theta);
            const double sin_t = std::sin(theta);
            const double A =
                -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
            const double B =
                -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
            // Translation to coordinate system R' where aligned paraboloid
            // valid Translation is R' = {x' = x + u, y' = y + v, z' = z + w}
            const double denominator = IRL::safelyTiny(4.0 * d * f - e * e);
            const double u = (2.0 * b * f - c * e) / denominator;
            const double v = -(b * e - 2.0 * d * c) / denominator;
            const double w =
                -(a + (-b * b * f + b * c * e - c * c * d) / denominator);

            IRL::UnitQuaternion rotation(theta, fit_frame[2]);
            IRL::Pt datum =
                pref - u * fit_frame[0] - v * fit_frame[1] - w * fit_frame[2];
            auto new_frame = rotation * fit_frame;
            const double max_curvature_dx = 1.0;
            double a_coeff = A;
            double b_coeff = B;
            if (std::sqrt(u * u + v * v + w * w) > 10.0 * mesh.dx() ||
                std::fabs(A) * mesh.dx() > max_curvature_dx ||
                std::fabs(B) * mesh.dx() > max_curvature_dx) {
              paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
            } else {
              paraboloid = IRL::Paraboloid(datum, new_frame, a_coeff, b_coeff);
            }
          } else {
            paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
          }

          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, a_liquid_volume_fraction(i, j, k), 1.0e-14,
                              paraboloid);
          if (solver_distance.getDistance() == -DBL_MAX) {
            (*a_interface)(i, j, k) =
                IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
          } else {
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            (*a_interface)(i, j, k) = paraboloid;
          }
        }
      }
    }
  }

  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void PLIC::getReconstruction(const Data<double>& a_liquid_volume_fraction,
                             const double a_dt, const Data<double>& a_U,
                             const Data<double>& a_V, const Data<double>& a_W,
                             Data<IRL::Paraboloid>* a_interface) {
  const BasicMesh& mesh = a_U.getMesh();

  Data<IRL::PlanarSeparator> interface(&mesh);
  updateReconstructionELVIRA(a_liquid_volume_fraction, &interface);
  updateReconstructionLVIRA(a_liquid_volume_fraction, 1, &interface);
  Data<IRL::Polygon> polygon(&mesh);
  updatePolygon(a_liquid_volume_fraction, interface, &polygon);
  polygon.updateBorder();

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] -= mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[0] += mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] -= mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[1] += mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] -= mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& pt : polygon(i, j, k)) {
          pt[2] += mesh.lz();
        }
      }
    }
  }

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else if (a_liquid_volume_fraction(i, j, k) >
                   IRL::global_constants::VF_HIGH) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          // continue;
        } else {
          const IRL::Normal norm_poly = polygon(i, j, k).calculateNormal();
          const IRL::Pt pref = polygon(i, j, k).calculateCentroid();
          IRL::ReferenceFrame fit_frame;
          int largest_dir = 0;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[1]))
            largest_dir = 1;
          if (std::fabs(norm_poly[largest_dir]) < std::fabs(norm_poly[2]))
            largest_dir = 2;
          if (largest_dir == 0)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 1)
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(0.0, 0.0, 1.0));
          else
            fit_frame[0] = crossProduct(norm_poly, IRL::Normal(1.0, 0.0, 0.0));
          fit_frame[0].normalize();
          fit_frame[1] = crossProduct(norm_poly, fit_frame[0]);
          fit_frame[2] = norm_poly;
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Paraboloid paraboloid;

          paraboloid = IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);

          auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                              upper_cell_pt);
          IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
              solver_distance(cell, a_liquid_volume_fraction(i, j, k), 1.0e-14,
                              paraboloid);

          if (solver_distance.getDistance() == -DBL_MAX) {
            (*a_interface)(i, j, k) =
                IRL::Paraboloid(pref, fit_frame, 1.0e-3, -1.0e-3);
          } else {
            auto new_datum =
                IRL::Pt(paraboloid.getDatum() +
                        solver_distance.getDistance() * fit_frame[2]);
            paraboloid.setDatum(new_datum);
            (*a_interface)(i, j, k) = paraboloid;
          }
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
