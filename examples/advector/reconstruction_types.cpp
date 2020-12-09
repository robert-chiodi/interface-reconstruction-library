// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/advector/reconstruction_types.h"

#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/interface_reconstruction_methods/elvira.h"
#include "irl/interface_reconstruction_methods/lvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/lvira_optimization.h"
#include "irl/interface_reconstruction_methods/r2p_neighborhood.h"
#include "irl/interface_reconstruction_methods/r2p_optimization.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/localizer_link_from_localized_separator_link.h"

#include "examples/advector/basic_mesh.h"
#include "examples/advector/data.h"
#include "examples/advector/vof_advection.h"

void getReconstruction(
    const std::string& a_reconstruction_method,
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedSeparatorLink>& a_localized_separator_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::PlanarSeparator>* a_interface) {
  if (a_reconstruction_method == "ELVIRA2D") {
    ELVIRA2D::getReconstruction(a_liquid_volume_fraction, a_dt, a_U, a_V, a_W,
                                a_interface);
  } else if (a_reconstruction_method == "LVIRA2D") {
    LVIRA2D::getReconstruction(a_liquid_volume_fraction, a_liquid_centroid,
                               a_gas_centroid, a_dt, a_U, a_V, a_W,
                               a_interface);
  } else if (a_reconstruction_method == "MOF2D") {
    MOF2D::getReconstruction(a_liquid_volume_fraction, a_liquid_centroid,
                             a_gas_centroid, a_localized_separator_link, a_dt,
                             a_U, a_V, a_W, a_interface);
  } else if (a_reconstruction_method == "AdvectedNormals") {
    AdvectedNormals::getReconstruction(
        a_liquid_volume_fraction, a_liquid_centroid, a_gas_centroid,
        a_localized_separator_link, a_dt, a_U, a_V, a_W, a_interface);
  } else if (a_reconstruction_method == "R2P2D") {
    R2P2D::getReconstruction(a_liquid_volume_fraction, a_liquid_centroid,
                             a_gas_centroid, a_localized_separator_link, a_dt,
                             a_U, a_V, a_W, a_interface);
  } else if (a_reconstruction_method == "ELVIRA3D") {
    ELVIRA3D::getReconstruction(a_liquid_volume_fraction, a_dt, a_U, a_V, a_W,
                                a_interface);
  } else if (a_reconstruction_method == "LVIRA3D") {
    LVIRA3D::getReconstruction(a_liquid_volume_fraction, a_liquid_centroid,
                               a_gas_centroid, a_dt, a_U, a_V, a_W,
                               a_interface);
  } else if (a_reconstruction_method == "MOF3D") {
    MOF3D::getReconstruction(a_liquid_volume_fraction, a_liquid_centroid,
                             a_gas_centroid, a_localized_separator_link, a_dt,
                             a_U, a_V, a_W, a_interface);
  } else if (a_reconstruction_method == "AdvectedNormals3D") {
    AdvectedNormals3D::getReconstruction(
        a_liquid_volume_fraction, a_liquid_centroid, a_gas_centroid,
        a_localized_separator_link, a_dt, a_U, a_V, a_W, a_interface);
  } else {
    std::cout << "Unknown reconstruction method of : "
              << a_reconstruction_method << '\n';
    std::cout << "Value entries are: ELVIRA2D, LVIRA2D, MOF2D, "
                 "AdvectedNormals, R2P2D, ELVIRA3D. \n";
    std::exit(-1);
  }
}

void ELVIRA2D::getReconstruction(const Data<double>& a_liquid_volume_fraction,
                                 const double a_dt, const Data<double>& a_U,
                                 const Data<double>& a_V,
                                 const Data<double>& a_W,
                                 Data<IRL::PlanarSeparator>* a_interface) {
  IRL::ELVIRANeighborhood neighborhood;
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  neighborhood.resize(9);
  IRL::RectangularCuboid cells[9];
  // Loop over cells in domain. Skip if cell is not mixed phase.
  const int k = 0;
  const int kk = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
          a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
        const double distance =
            std::copysign(IRL::global_constants::ARBITRARILY_LARGE_DISTANCE,
                          a_liquid_volume_fraction(i, j, k) - 0.5);
        (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
            IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
        continue;
      }
      // Build surrounding stencil information for ELVIRA.
      for (int ii = i - 1; ii < i + 2; ++ii) {
        for (int jj = j - 1; jj < j + 2; ++jj) {
          // Reversed order, bad for cache locality but thats okay..
          cells[(jj - j + 1) * 3 + (ii - i + 1)] =
              IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
          neighborhood.setMember(&cells[(jj - j + 1) * 3 + (ii - i + 1)],
                                 &a_liquid_volume_fraction(ii, jj, 0), ii - i,
                                 jj - j);
        }
      }
      // Now perform actual ELVIRA and obtain interface PlanarSeparator
      (*a_interface)(i, j, k) = reconstructionWithELVIRA2D(neighborhood);
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void ELVIRA3D::getReconstruction(const Data<double>& a_liquid_volume_fraction,
                                 const double a_dt, const Data<double>& a_U,
                                 const Data<double>& a_V,
                                 const Data<double>& a_W,
                                 Data<IRL::PlanarSeparator>* a_interface) {
  IRL::ELVIRANeighborhood neighborhood;
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  neighborhood.resize(27);
  IRL::RectangularCuboid cells[27];
  // Loop over cells in domain. Skip if cell is not mixed phase.
  // const int k = 0;
  // const int kk = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
            a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, a_liquid_volume_fraction(i, j, k) - 0.5);
          (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }
        // Build surrounding stencil information for ELVIRA.
        for (int ii = i - 1; ii < i + 2; ++ii) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            for (int kk = k - 1; kk < k + 2; ++kk) {
              // Reversed order, bad for cache locality but thats okay..
              cells[9 * (kk - k + 1 ) + (jj - j + 1) * 3 + (ii - i + 1)] =
                  IRL::RectangularCuboid::fromBoundingPts(
                      IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                      IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              neighborhood.setMember(&cells[9 * (kk - k + 1 ) + (jj - j + 1) * 3 + (ii - i + 1)],
                                    &a_liquid_volume_fraction(ii, jj, kk), ii - i,
                                    jj - j, kk - k);
            }
          }
        }
        // Now perform actual ELVIRA and obtain interface PlanarSeparator
        (*a_interface)(i, j, k) = reconstructionWithELVIRA3D(neighborhood);
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void LVIRA2D::getReconstruction(const Data<double>& a_liquid_volume_fraction,
                                const Data<IRL::Pt>& a_liquid_centroid,
                                const Data<IRL::Pt>& a_gas_centroid,
                                const double a_dt, const Data<double>& a_U,
                                const Data<double>& a_V,
                                const Data<double>& a_W,
                                Data<IRL::PlanarSeparator>* a_interface) {
  IRL::LVIRANeighborhood<IRL::RectangularCuboid> neighborhood;
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  neighborhood.resize(9);
  neighborhood.setCenterOfStencil(4);
  IRL::RectangularCuboid cells[9];
  // Loop over cells in domain. Skip if cell is not mixed phase.
  const int k = 0;
  const int kk = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
          a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
        const double distance =
            std::copysign(IRL::global_constants::ARBITRARILY_LARGE_DISTANCE,
                          a_liquid_volume_fraction(i, j, k) - 0.5);
        (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
            IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
        continue;
      }
      // Build surrounding stencil information for ELVIRA.
      for (int ii = i - 1; ii < i + 2; ++ii) {
        for (int jj = j - 1; jj < j + 2; ++jj) {
          // Reversed order, bad for cache locality but thats okay..
          cells[(jj - j + 1) * 3 + (ii - i + 1)] =
              IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
          neighborhood.setMember(static_cast<IRL::UnsignedIndex_t>(
                                     (jj - j + 1) * 3 + (ii - i + 1)),
                                 &cells[(jj - j + 1) * 3 + (ii - i + 1)],
                                 &a_liquid_volume_fraction(ii, jj, kk));
        }
      }
      // Now create initial guess using centroids
      auto bary_normal = IRL::Normal::fromPtNormalized(
          a_gas_centroid(i, j, k) - a_liquid_centroid(i, j, k));
      bary_normal[2] = 0.0;
      bary_normal.normalize();
      const double initial_distance =
          bary_normal * neighborhood.getCenterCell().calculateCentroid();
      (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
          IRL::Plane(bary_normal, initial_distance));
      setDistanceToMatchVolumeFractionPartialFill(
          neighborhood.getCenterCell(),
          neighborhood.getCenterCellStoredMoments(), &(*a_interface)(i, j, k));

      (*a_interface)(i, j, k) =
          reconstructionWithLVIRA2D(neighborhood, (*a_interface)(i, j, k));
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void LVIRA3D::getReconstruction(const Data<double>& a_liquid_volume_fraction,
                                const Data<IRL::Pt>& a_liquid_centroid,
                                const Data<IRL::Pt>& a_gas_centroid,
                                const double a_dt, const Data<double>& a_U,
                                const Data<double>& a_V,
                                const Data<double>& a_W,
                                Data<IRL::PlanarSeparator>* a_interface) {
  IRL::LVIRANeighborhood<IRL::RectangularCuboid> neighborhood;
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  neighborhood.resize(27);
  neighborhood.setCenterOfStencil(13);
  IRL::RectangularCuboid cells[27];
  // Loop over cells in domain. Skip if cell is not mixed phase.
  // const int k = 0;
  // const int kk = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
            a_liquid_volume_fraction(i, j, k) >
                IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, a_liquid_volume_fraction(i, j, k) - 0.5);
          (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }
        // Build surrounding stencil information for ELVIRA.
        for (int ii = i - 1; ii < i + 2; ++ii) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            for (int kk = k - 1; kk < k + 2; ++kk) {
              // Reversed order, bad for cache locality but thats okay..
              const int local_index =
                  (kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1);
              cells[local_index] = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                  IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
              neighborhood.setMember(
                  static_cast<IRL::UnsignedIndex_t>(local_index),
                  &cells[local_index], &a_liquid_volume_fraction(ii, jj, kk));
            }
          }
        }
        // Now create initial guess using centroids
        auto bary_normal = IRL::Normal::fromPtNormalized(
            a_gas_centroid(i, j, k) - a_liquid_centroid(i, j, k));
        bary_normal.normalize();
        const double initial_distance =
            bary_normal * neighborhood.getCenterCell().calculateCentroid();
        (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
            IRL::Plane(bary_normal, initial_distance));
        setDistanceToMatchVolumeFractionPartialFill(
            neighborhood.getCenterCell(),
            neighborhood.getCenterCellStoredMoments(),
            &(*a_interface)(i, j, k));

        (*a_interface)(i, j, k) =
            reconstructionWithLVIRA3D(neighborhood, (*a_interface)(i, j, k));
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void MOF2D::getReconstruction(
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedSeparatorLink>& a_localized_separator_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::PlanarSeparator>* a_interface) {
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();

  const int k = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
          a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
        const double distance =
            std::copysign(IRL::global_constants::ARBITRARILY_LARGE_DISTANCE,
                          a_liquid_volume_fraction(i, j, k) - 0.5);
        (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
            IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
        continue;
      }
      auto cell = IRL::RectangularCuboid::fromBoundingPts(
          IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
          IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
      double vol = cell.calculateVolume();
      IRL::SeparatedMoments<IRL::VolumeMoments> svm(
          IRL::VolumeMoments(a_liquid_volume_fraction(i, j, k) * vol,
                             a_liquid_centroid(i, j, k)),
          IRL::VolumeMoments((1.0 - a_liquid_volume_fraction(i, j, k)) * vol,
                             a_gas_centroid(i, j, k)));
      (*a_interface)(i, j, k) =
          IRL::reconstructionWithMOF2D(cell, svm, 0.5, 0.5);
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void MOF3D::getReconstruction(
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedSeparatorLink>& a_localized_separator_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::PlanarSeparator>* a_interface) {
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();

  // const int k = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
            a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, a_liquid_volume_fraction(i, j, k) - 0.5);
          (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        }
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        double vol = cell.calculateVolume();
        IRL::SeparatedMoments<IRL::VolumeMoments> svm(
            IRL::VolumeMoments(a_liquid_volume_fraction(i, j, k) * vol,
                              a_liquid_centroid(i, j, k)),
            IRL::VolumeMoments((1.0 - a_liquid_volume_fraction(i, j, k)) * vol,
                              a_gas_centroid(i, j, k)));
        (*a_interface)(i, j, k) =
            IRL::reconstructionWithMOF3D(cell, svm, 0.5, 0.5);
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void AdvectedNormals::getReconstruction(
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedSeparatorLink>& a_localized_separator_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::PlanarSeparator>* a_interface) {
  // Get mesh everything is living on.
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  // Container for moments from advection
  Data<IRL::ListedVolumeMoments<IRL::VolumeMomentsAndNormal>> listed_moments(
      &mesh);

  const int k = 0;
  const int kk = 0;
  for (int i = mesh.imino() + 1; i <= mesh.imaxo() - 1; ++i) {
    for (int j = mesh.jmino() + 1; j <= mesh.jmaxo() - 1; ++j) {
      auto cell = IRL::RectangularCuboid::fromBoundingPts(
          IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
          IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
      const auto localizer_link = IRL::LocalizerLinkFromLocalizedSeparatorLink(
          &a_localized_separator_link(i, j, k));
      for (IRL::UnsignedIndex_t n = 0;
           n < (*a_interface)(i, j, k).getNumberOfPlanes(); ++n) {
        IRL::Polygon interface_poly =
            IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                cell, (*a_interface)(i, j, k), (*a_interface)(i, j, k)[n]);
        if (interface_poly.getNumberOfVertices() == 0) {
          continue;
        }
        for (IRL::UnsignedIndex_t tri = 0;
             tri < interface_poly.getNumberOfSimplicesInDecomposition();
             ++tri) {
          IRL::Tri simplex = static_cast<IRL::Tri>(
              interface_poly.getSimplexFromDecomposition(tri));
          for (auto& vertex : simplex) {
            vertex = back_project_vertex(vertex, a_dt, a_U, a_V, a_W);
          }
          simplex.calculateAndSetPlaneOfExistence();
          auto new_moments =
              IRL::getVolumeMoments<IRL::TaggedAccumulatedListedVolumeMoments<
                  IRL::VolumeMomentsAndNormal>>(simplex, localizer_link);
          for (IRL::UnsignedIndex_t moment = 0; moment < new_moments.size();
               ++moment) {
            auto index_for_tag =
                getIndexFromTag(mesh, new_moments.getTagForIndex(moment));
            listed_moments(index_for_tag[0], index_for_tag[1],
                           index_for_tag[2]) +=
                new_moments.getMomentsForIndex(moment);
          }
        }
      }
    }
  }

  // Remove Z components from advected surface elements that will be used
  // This can occur by polygons being rotated by the flow field during
  // advection.
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      const IRL::UnsignedIndex_t starting_length =
          listed_moments(i, j, k).size();
      for (IRL::UnsignedIndex_t n = starting_length - 1;
           n != static_cast<IRL::UnsignedIndex_t>(-1); --n) {
        IRL::VolumeMomentsAndNormal& moment = listed_moments(i, j, k)[n];
        moment.normalizeByVolume();
        moment.normal()[2] = 0.0;
        moment.normal().normalize();
        if (moment.normal().calculateMagnitude() < 0.95) {
          listed_moments(i, j, k).erase(n);
        } else {
          moment.multiplyByVolume();
        }
      }
    }
  }

  // Now have all of the advected moments. Get and store the reconstructions
  // from this by partitioning with Kmeans.
  IRL::R2PNeighborhood<IRL::RectangularCuboid> neighborhood;
  ////////////////////////////////////////////////
  neighborhood.resize(9);  // SET TO 9 BECAUSE OF 2D.
  neighborhood.setCenterOfStencil(4);
  IRL::RectangularCuboid stencil_cells[9];
  IRL::SeparatedMoments<IRL::VolumeMoments> stencil_moments[9];
  int num_mof = 0;
  int num_adv = 0;
  int num_adv2 = 0;
  ////////////////////////////////////////////////
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
          a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
        const double distance =
            std::copysign(IRL::global_constants::ARBITRARILY_LARGE_DISTANCE,
                          a_liquid_volume_fraction(i, j, k) - 0.5);
        (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
            IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
        continue;
      } else if (listed_moments(i, j, k).size() == 0) {
        // No interface advected in, use MoF
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        double vol = cell.calculateVolume();
        IRL::SeparatedMoments<IRL::VolumeMoments> svm(
            IRL::VolumeMoments(a_liquid_volume_fraction(i, j, k) * vol,
                               a_liquid_centroid(i, j, k)),
            IRL::VolumeMoments((1.0 - a_liquid_volume_fraction(i, j, k)) * vol,
                               a_gas_centroid(i, j, k)));
        (*a_interface)(i, j, k) = IRL::reconstructionWithMOF2D(cell, svm);
        ++num_mof;
      } else {
        // Set up R2P neighborhood
        for (int ii = i - 1; ii < i + 2; ++ii) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            const int ind = (ii - i + 1) * 3 + (jj - j + 1);
            stencil_cells[ind] = IRL::RectangularCuboid::fromBoundingPts(
                IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(k)),
                IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(k + 1)));
            double vol = stencil_cells[ind].calculateVolume();
            stencil_moments[ind] = IRL::SeparatedMoments<IRL::VolumeMoments>(
                IRL::VolumeMoments(a_liquid_volume_fraction(ii, jj, kk) * vol,
                                   a_liquid_centroid(ii, jj, kk)),
                IRL::VolumeMoments(
                    (1.0 - a_liquid_volume_fraction(ii, jj, kk)) * vol,
                    a_gas_centroid(ii, jj, kk)));
            neighborhood.setMember(static_cast<IRL::UnsignedIndex_t>(ind),
                                   &stencil_cells[ind], &stencil_moments[ind]);
          }
        }
        (*a_interface)(i, j, k) = IRL::reconstructionWithAdvectedNormals(
            listed_moments(i, j, k), neighborhood);
        ++num_adv;
        if ((*a_interface)(i, j, k).getNumberOfPlanes() == 2) {
          ++num_adv2;
        }
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void AdvectedNormals3D::getReconstruction(
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedSeparatorLink>& a_localized_separator_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::PlanarSeparator>* a_interface) {
  // Get mesh everything is living on.
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  // Container for moments from advection
  Data<IRL::ListedVolumeMoments<IRL::VolumeMomentsAndNormal>> listed_moments(
      &mesh);

  // const int k = 0;
  // const int kk = 0;
  for (int i = mesh.imino() + 1; i <= mesh.imaxo() - 1; ++i) {
    for (int j = mesh.jmino() + 1; j <= mesh.jmaxo() - 1; ++j) {
      for (int k = mesh.kmino() + 1; k <= mesh.kmaxo() - 1; ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto localizer_link = IRL::LocalizerLinkFromLocalizedSeparatorLink(
            &a_localized_separator_link(i, j, k));
        for (IRL::UnsignedIndex_t n = 0;
            n < (*a_interface)(i, j, k).getNumberOfPlanes(); ++n) {
          IRL::Polygon interface_poly =
              IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                  cell, (*a_interface)(i, j, k), (*a_interface)(i, j, k)[n]);
          if (interface_poly.getNumberOfVertices() == 0) {
            continue;
          }
          for (IRL::UnsignedIndex_t tri = 0;
              tri < interface_poly.getNumberOfSimplicesInDecomposition();
              ++tri) {
            IRL::Tri simplex = static_cast<IRL::Tri>(
                interface_poly.getSimplexFromDecomposition(tri));
            for (auto& vertex : simplex) {
              vertex = back_project_vertex(vertex, a_dt, a_U, a_V, a_W);
            }
            simplex.calculateAndSetPlaneOfExistence();
            auto new_moments =
                IRL::getVolumeMoments<IRL::TaggedAccumulatedListedVolumeMoments<
                    IRL::VolumeMomentsAndNormal>>(simplex, localizer_link);
            for (IRL::UnsignedIndex_t moment = 0; moment < new_moments.size();
                ++moment) {
              auto index_for_tag =
                  getIndexFromTag(mesh, new_moments.getTagForIndex(moment));
              listed_moments(index_for_tag[0], index_for_tag[1],
                            index_for_tag[2]) +=
                  new_moments.getMomentsForIndex(moment);
            }
          }
        }
      }
    }
  }

  // Remove Z components from advected surface elements that will be used
  // This can occur by polygons being rotated by the flow field during
  // advection.
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        const IRL::UnsignedIndex_t starting_length =
            listed_moments(i, j, k).size();
        for (IRL::UnsignedIndex_t n = starting_length - 1;
            n != static_cast<IRL::UnsignedIndex_t>(-1); --n) {
          IRL::VolumeMomentsAndNormal& moment = listed_moments(i, j, k)[n];
          moment.normalizeByVolume();
          // moment.normal()[2] = 0.0;
          moment.normal().normalize();
          if (moment.normal().calculateMagnitude() < 0.95) {
            listed_moments(i, j, k).erase(n);
          } else {
            moment.multiplyByVolume();
          }
        }
      }
    }
  }

  // Now have all of the advected moments. Get and store the reconstructions
  // from this by partitioning with Kmeans.
  IRL::R2PNeighborhood<IRL::RectangularCuboid> neighborhood;
  ////////////////////////////////////////////////
  neighborhood.resize(27);  // SET TO 9 BECAUSE OF 2D.
  neighborhood.setCenterOfStencil(13);
  IRL::RectangularCuboid stencil_cells[27];
  IRL::SeparatedMoments<IRL::VolumeMoments> stencil_moments[27];
  int num_mof = 0;
  int num_adv = 0;
  int num_adv2 = 0;
  ////////////////////////////////////////////////
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
            a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
          const double distance =
              std::copysign(1.0, a_liquid_volume_fraction(i, j, k) - 0.5);
          (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          continue;
        } else if (listed_moments(i, j, k).size() == 0) {
          // No interface advected in, use MoF
          auto cell = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
              IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
          double vol = cell.calculateVolume();
          IRL::SeparatedMoments<IRL::VolumeMoments> svm(
              IRL::VolumeMoments(a_liquid_volume_fraction(i, j, k) * vol,
                                a_liquid_centroid(i, j, k)),
              IRL::VolumeMoments((1.0 - a_liquid_volume_fraction(i, j, k)) * vol,
                                a_gas_centroid(i, j, k)));
          (*a_interface)(i, j, k) = IRL::reconstructionWithMOF3D(cell, svm);
          ++num_mof;
        } else {
          // Set up R2P neighborhood
          for (int ii = i - 1; ii < i + 2; ++ii) {
            for (int jj = j - 1; jj < j + 2; ++jj) {
              for (int kk = k - 1; kk < k + 2; ++kk) {
                const int ind = (ii - i + 1) * 9 + (jj - j + 1) * 3 + (kk - k + 1);
                stencil_cells[ind] = IRL::RectangularCuboid::fromBoundingPts(
                    IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(kk)),
                    IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(kk + 1)));
                double vol = stencil_cells[ind].calculateVolume();
                stencil_moments[ind] = IRL::SeparatedMoments<IRL::VolumeMoments>(
                    IRL::VolumeMoments(a_liquid_volume_fraction(ii, jj, kk) * vol,
                                      a_liquid_centroid(ii, jj, kk)),
                    IRL::VolumeMoments(
                        (1.0 - a_liquid_volume_fraction(ii, jj, kk)) * vol,
                        a_gas_centroid(ii, jj, kk)));
                neighborhood.setMember(static_cast<IRL::UnsignedIndex_t>(ind),
                                      &stencil_cells[ind], &stencil_moments[ind]);
              }
            }
          }
          (*a_interface)(i, j, k) = IRL::reconstructionWithAdvectedNormals(
              listed_moments(i, j, k), neighborhood);
          ++num_adv;
          if ((*a_interface)(i, j, k).getNumberOfPlanes() == 2) {
            ++num_adv2;
          }
        }
      }
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void R2P2D::getReconstruction(
    const Data<double>& a_liquid_volume_fraction,
    const Data<IRL::Pt>& a_liquid_centroid, const Data<IRL::Pt>& a_gas_centroid,
    const Data<IRL::LocalizedSeparatorLink>& a_localized_separator_link,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::PlanarSeparator>* a_interface) {
  // Get mesh everything is living on.
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  // Container for moments from advection
  Data<IRL::ListedVolumeMoments<IRL::VolumeMomentsAndNormal>> listed_moments(
      &mesh);

  const int k = 0;
  const int kk = 0;
  for (int i = mesh.imino() + 1; i <= mesh.imaxo() - 1; ++i) {
    for (int j = mesh.jmino() + 1; j <= mesh.jmaxo() - 1; ++j) {
      auto cell = IRL::RectangularCuboid::fromBoundingPts(
          IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
          IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
      const auto localizer_link = IRL::LocalizerLinkFromLocalizedSeparatorLink(
          &a_localized_separator_link(i, j, k));
      for (IRL::UnsignedIndex_t n = 0;
           n < (*a_interface)(i, j, k).getNumberOfPlanes(); ++n) {
        IRL::Polygon interface_poly =
            IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                cell, (*a_interface)(i, j, k), (*a_interface)(i, j, k)[n]);
        if (interface_poly.getNumberOfVertices() == 0) {
          continue;
        }
        for (IRL::UnsignedIndex_t tri = 0;
             tri < interface_poly.getNumberOfSimplicesInDecomposition();
             ++tri) {
          IRL::Tri simplex = static_cast<IRL::Tri>(
              interface_poly.getSimplexFromDecomposition(tri));
          for (auto& vertex : simplex) {
            vertex = back_project_vertex(vertex, a_dt, a_U, a_V, a_W);
          }
          simplex.calculateAndSetPlaneOfExistence();
          auto new_moments =
              IRL::getVolumeMoments<IRL::TaggedAccumulatedListedVolumeMoments<
                  IRL::VolumeMomentsAndNormal>>(simplex, localizer_link);
          for (IRL::UnsignedIndex_t moment = 0; moment < new_moments.size();
               ++moment) {
            auto index_for_tag =
                getIndexFromTag(mesh, new_moments.getTagForIndex(moment));
            listed_moments(index_for_tag[0], index_for_tag[1],
                           index_for_tag[2]) +=
                new_moments.getMomentsForIndex(moment);
          }
        }
      }
    }
  }

  // Remove Z components from advected surface elements that will be used
  // This can occur by polygons being rotated by the flow field during
  // advection.
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      const IRL::UnsignedIndex_t starting_length =
          listed_moments(i, j, k).size();
      for (IRL::UnsignedIndex_t n = starting_length - 1;
           n != static_cast<IRL::UnsignedIndex_t>(-1); --n) {
        IRL::VolumeMomentsAndNormal& moment = listed_moments(i, j, k)[n];
        moment.normalizeByVolume();
        moment.normal()[2] = 0.0;
        moment.normal().normalize();
        if (moment.normal().calculateMagnitude() < 0.95) {
          listed_moments(i, j, k).erase(n);
        } else {
          moment.multiplyByVolume();
        }
      }
    }
  }

  // Now have all of the advected moments. Get and store the reconstructions
  // from this by partitioning with Kmeans.
  IRL::R2PNeighborhood<IRL::RectangularCuboid> neighborhood;
  ////////////////////////////////////////////////
  neighborhood.resize(9);  // SET TO 9 BECAUSE OF 2D.
  neighborhood.setCenterOfStencil(4);
  IRL::RectangularCuboid stencil_cells[9];
  IRL::SeparatedMoments<IRL::VolumeMoments> stencil_moments[9];
  int num_mof = 0;
  int num_adv = 0;
  int num_adv2 = 0;
  ////////////////////////////////////////////////
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (a_liquid_volume_fraction(i, j, k) < IRL::global_constants::VF_LOW ||
          a_liquid_volume_fraction(i, j, k) > IRL::global_constants::VF_HIGH) {
        const double distance =
            std::copysign(IRL::global_constants::ARBITRARILY_LARGE_DISTANCE,
                          a_liquid_volume_fraction(i, j, k) - 0.5);
        (*a_interface)(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
            IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
        continue;
      }

      // Set up R2P neighborhood
      for (int ii = i - 1; ii < i + 2; ++ii) {
        for (int jj = j - 1; jj < j + 2; ++jj) {
          const int ind = (ii - i + 1) * 3 + (jj - j + 1);
          stencil_cells[ind] = IRL::RectangularCuboid::fromBoundingPts(
              IRL::Pt(mesh.x(ii), mesh.y(jj), mesh.z(k)),
              IRL::Pt(mesh.x(ii + 1), mesh.y(jj + 1), mesh.z(k + 1)));
          double vol = stencil_cells[ind].calculateVolume();
          stencil_moments[ind] = IRL::SeparatedMoments<IRL::VolumeMoments>(
              IRL::VolumeMoments(a_liquid_volume_fraction(ii, jj, kk) * vol,
                                 a_liquid_centroid(ii, jj, kk)),
              IRL::VolumeMoments(
                  (1.0 - a_liquid_volume_fraction(ii, jj, kk)) * vol,
                  a_gas_centroid(ii, jj, kk)));
          neighborhood.setMember(static_cast<IRL::UnsignedIndex_t>(ind),
                                 &stencil_cells[ind], &stencil_moments[ind]);
        }
      }

      if (listed_moments(i, j, k).size() == 0) {
        // No interface advected in, use MoF
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        double vol = cell.calculateVolume();
        IRL::SeparatedMoments<IRL::VolumeMoments> svm(
            IRL::VolumeMoments(a_liquid_volume_fraction(i, j, k) * vol,
                               a_liquid_centroid(i, j, k)),
            IRL::VolumeMoments((1.0 - a_liquid_volume_fraction(i, j, k)) * vol,
                               a_gas_centroid(i, j, k)));
        (*a_interface)(i, j, k) = IRL::reconstructionWithMOF2D(cell, svm);
        ++num_mof;
        neighborhood.setSurfaceArea(
            getReconstructionSurfaceArea(cell, (*a_interface)(i, j, k)));
      } else {
        (*a_interface)(i, j, k) = IRL::reconstructionWithAdvectedNormals(
            listed_moments(i, j, k), neighborhood);
        ++num_adv;
        if ((*a_interface)(i, j, k).getNumberOfPlanes() == 2) {
          ++num_adv2;
        }
        double area_sum = 0.0;
        for (const auto& moment : listed_moments(i, j, k)) {
          area_sum += moment.volumeMoments().volume();
        }
        neighborhood.setSurfaceArea(area_sum);
      }
      (*a_interface)(i, j, k) =
          reconstructionWithR2P2D(neighborhood, (*a_interface)(i, j, k));
    }
  }
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void correctInterfacePlaneBorders(Data<IRL::PlanarSeparator>* a_interface) {
  const BasicMesh& mesh = (*a_interface).getMesh();
  // Fix distance to recreate volume fraction

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& plane : (*a_interface)(i, j, k)) {
          plane.distance() = plane.distance() - plane.normal()[0] * mesh.lx();
        }
      }
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& plane : (*a_interface)(i, j, k)) {
          plane.distance() = plane.distance() + plane.normal()[0] * mesh.lx();
        }
      }
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& plane : (*a_interface)(i, j, k)) {
          plane.distance() = plane.distance() - plane.normal()[1] * mesh.ly();
        }
      }
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        for (auto& plane : (*a_interface)(i, j, k)) {
          plane.distance() = plane.distance() + plane.normal()[1] * mesh.ly();
        }
      }
    }
  }

  // z- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k < mesh.kmin(); ++k) {
        for (auto& plane : (*a_interface)(i, j, k)) {
          plane.distance() = plane.distance() - plane.normal()[2] * mesh.lz();
        }
      }
    }
  }

  // z+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmax() + 1; k <= mesh.kmaxo(); ++k) {
        for (auto& plane : (*a_interface)(i, j, k)) {
          plane.distance() = plane.distance() - plane.normal()[2] * mesh.lz();
        }
      }
    }
  }
}
