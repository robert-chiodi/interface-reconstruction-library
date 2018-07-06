// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/interface_reconstruction_methods/elvira.h"

#include <cmath>
#include <iostream>

#include "src/geometry/polygons/polygon.h"
#include "src/interface_reconstruction_methods/plane_distance.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

PlanarSeparator ELVIRA_2D::solve(
    const ELVIRANeighborhood* a_neighborhood_pointer) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_VF_m = a_neighborhood_pointer;
  return this->solve();
}

PlanarSeparator ELVIRA_2D::solve(void) {
  minimum_error_m = DBL_MAX;
  guess_reconstruction_m =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));
  double tmp_normal_approximations[3];

  // Assume x integration direction
  this->makeColumnSumsX();
  double assumed_normal = neighborhood_VF_m->getStoredMoments(1, 0) >
                                  neighborhood_VF_m->getStoredMoments(-1, 0)
                              ? -1.0
                              : 1.0;
  this->fillOtherNormals(1, tmp_normal_approximations);
  // Try all normal combinations
  for (int y_normal = 0; y_normal < 3; ++y_normal) {
    this->tryNormal(
        Normal(assumed_normal, tmp_normal_approximations[y_normal], 0.0));
  }

  // Assume y integration direction
  this->makeColumnSumsY();
  assumed_normal = neighborhood_VF_m->getStoredMoments(0, 1) >
                           neighborhood_VF_m->getStoredMoments(0, -1)
                       ? -1.0
                       : 1.0;
  this->fillOtherNormals(0, tmp_normal_approximations);

  // Try all normal combinations
  for (int x_normal = 0; x_normal < 3; ++x_normal) {
    this->tryNormal(
        Normal(tmp_normal_approximations[x_normal], assumed_normal, 0.0));
  }

  return best_reconstruction_m;
}

double& ELVIRA_2D::columnSum(const int i) {
  return column_sums_m[calculateLinearIndex(i)];
}

Pt& ELVIRA_2D::columnCenters(const int i) {
  return column_centers_m[calculateLinearIndex(i)];
}

double ELVIRA_2D::columnCenters(const UnsignedIndex_t dimension, const int i) {
  assert(dimension < 2);
  return column_centers_m[calculateLinearIndex(i)][dimension];
}

UnsignedIndex_t ELVIRA_2D::calculateLinearIndex(const int i) {
  assert(i >= -1 && i < 2);
  return static_cast<UnsignedIndex_t>(i + 1);
}

void ELVIRA_2D::fillOtherNormals(const UnsignedIndex_t a_dimension_index,
                                 double a_normal_approximation_holder[3]) {
  // First index: backward, central, forward
  a_normal_approximation_holder[0] =
      this->computeDerivative0(0, -1, a_dimension_index);
  a_normal_approximation_holder[1] =
      this->computeDerivative0(1, -1, a_dimension_index);
  a_normal_approximation_holder[2] =
      this->computeDerivative0(1, 0, a_dimension_index);
}

double ELVIRA_2D::computeDerivative0(const int a_left_index,
                                     const int a_right_index,
                                     const UnsignedIndex_t spatial_index) {
  return -(this->columnSum(a_left_index) - this->columnSum(a_right_index)) /
         (columnCenters(spatial_index, a_left_index) -
          columnCenters(spatial_index, a_right_index));
}

void ELVIRA_2D::makeColumnSumsX(void) {
  for (int j = -1; j < 2; ++j) {
    double sum = {0.0};
    for (int i = -1; i < 2; ++i) {
      sum += neighborhood_VF_m->getStoredMoments(i, j) *
             (neighborhood_VF_m->getCell(i, j)).calculateSideLength(0);
    }
    columnCenters(j) = (neighborhood_VF_m->getCell(0, j)).calculateCentroid();
    this->columnSum(j) = sum;
  }
}

void ELVIRA_2D::makeColumnSumsY(void) {
  for (int i = -1; i < 2; ++i) {
    double sum = {0.0};
    for (int j = -1; j < 2; ++j) {
      sum += neighborhood_VF_m->getStoredMoments(i, j) *
             (neighborhood_VF_m->getCell(i, j)).calculateSideLength(1);
    }
    columnCenters(i) = (neighborhood_VF_m->getCell(i, 0)).calculateCentroid();
    this->columnSum(i) = sum;
  }
}

void ELVIRA_2D::tryNormal(Normal a_normal) {
  a_normal.normalize();
  guess_reconstruction_m[0] =
      Plane(a_normal, findDistanceOnePlane(
                          neighborhood_VF_m->getCell(0, 0),
                          neighborhood_VF_m->getStoredMoments(0, 0), a_normal));
  double error = {0.0};
  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      double volume_fraction =
          getVolumeFraction<ReconstructionDefaultCuttingMethod>(
              neighborhood_VF_m->getCell(i, j), guess_reconstruction_m);
      error += std::pow(
          (volume_fraction - neighborhood_VF_m->getStoredMoments(i, j)), 2);
    }
  }
  if (error < minimum_error_m) {
    best_reconstruction_m = guess_reconstruction_m;
    minimum_error_m = error;
  }
}

//******************************************************************* //
//     Function implementations for ELVIRA_3D below this.
//******************************************************************* //

PlanarSeparator ELVIRA_3D::solve(
    const ELVIRANeighborhood* a_neighborhood_pointer) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_VF_m = a_neighborhood_pointer;
  return this->solve();
}

PlanarSeparator ELVIRA_3D::solve(void) {
  minimum_error_m = DBL_MAX;
  guess_reconstruction_m =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));
  double tmp_normal_approximations[3][2];

  // Assume x integration direction
  this->makeColumnSumsX();
  double assumed_normal = neighborhood_VF_m->getStoredMoments(1, 0, 0) >
                                  neighborhood_VF_m->getStoredMoments(-1, 0, 0)
                              ? -1.0
                              : 1.0;
  this->fillOtherNormals(1, 2, tmp_normal_approximations);
  // Try all normal combinations
  for (int z_normal = 0; z_normal < 3; ++z_normal) {
    for (int y_normal = 0; y_normal < 3; ++y_normal) {
      this->tryNormal(Normal(assumed_normal,
                             tmp_normal_approximations[y_normal][0],
                             tmp_normal_approximations[z_normal][1]));
    }
  }

  // Assume y integration direction
  this->makeColumnSumsY();
  assumed_normal = neighborhood_VF_m->getStoredMoments(0, 1, 0) >
                           neighborhood_VF_m->getStoredMoments(0, -1, 0)
                       ? -1.0
                       : 1.0;
  this->fillOtherNormals(0, 2, tmp_normal_approximations);

  // Try all normal combinations
  for (int z_normal = 0; z_normal < 3; ++z_normal) {
    for (int x_normal = 0; x_normal < 3; ++x_normal) {
      this->tryNormal(Normal(tmp_normal_approximations[x_normal][0],
                             assumed_normal,
                             tmp_normal_approximations[z_normal][1]));
    }
  }

  // Assume z integration direction
  this->makeColumnSumsZ();
  assumed_normal = neighborhood_VF_m->getStoredMoments(0, 0, 1) >
                           neighborhood_VF_m->getStoredMoments(0, 0, -1)
                       ? -1.0
                       : 1.0;
  this->fillOtherNormals(0, 1, tmp_normal_approximations);

  // Try all normal combinations
  for (int y_normal = 0; y_normal < 3; ++y_normal) {
    for (int x_normal = 0; x_normal < 3; ++x_normal) {
      this->tryNormal(Normal(tmp_normal_approximations[x_normal][0],
                             tmp_normal_approximations[y_normal][1],
                             assumed_normal));
    }
  }
  return best_reconstruction_m;
}

double& ELVIRA_3D::columnSum(const int i, const int j) {
  return column_sums_m[calculateLinearIndex(i, j)];
}

Pt& ELVIRA_3D::columnCenters(const int i, const int j) {
  return column_centers_m[calculateLinearIndex(i, j)];
}

double ELVIRA_3D::columnCenters(const UnsignedIndex_t dimension, const int i,
                                const int j) {
  assert(dimension < 3);
  return (column_centers_m[calculateLinearIndex(i, j)])[dimension];
}

UnsignedIndex_t ELVIRA_3D::calculateLinearIndex(const int i, const int j) {
  assert(i >= -1 && i < 2);
  assert(j >= -1 && j < 2);
  return static_cast<UnsignedIndex_t>((i + 1) + (j + 1) * 3);
}

void ELVIRA_3D::fillOtherNormals(const UnsignedIndex_t a_first_column_index,
                                 const UnsignedIndex_t a_second_column_index,
                                 double a_normal_approximation_holder[3][2]) {
  // First index: backward, central, forward
  a_normal_approximation_holder[0][0] =
      this->computeDerivative0(0, -1, a_first_column_index);
  a_normal_approximation_holder[1][0] =
      this->computeDerivative0(1, -1, a_first_column_index);
  a_normal_approximation_holder[2][0] =
      this->computeDerivative0(1, 0, a_first_column_index);

  // Second index: backward, central, forward
  a_normal_approximation_holder[0][1] =
      this->computeDerivative1(0, -1, a_second_column_index);
  a_normal_approximation_holder[1][1] =
      this->computeDerivative1(1, -1, a_second_column_index);
  a_normal_approximation_holder[2][1] =
      this->computeDerivative1(1, 0, a_second_column_index);
}

double ELVIRA_3D::computeDerivative0(const int a_left_index,
                                     const int a_right_index,
                                     const UnsignedIndex_t a_spatial_index) {
  return -(this->columnSum(a_left_index, 0) -
           this->columnSum(a_right_index, 0)) /
         (columnCenters(a_spatial_index, a_left_index, 0) -
          columnCenters(a_spatial_index, a_right_index, 0));
}

double ELVIRA_3D::computeDerivative1(const int a_left_index,
                                     const int a_right_index,
                                     const UnsignedIndex_t a_spatial_index) {
  return -(this->columnSum(0, a_left_index) -
           this->columnSum(0, a_right_index)) /
         (columnCenters(a_spatial_index, 0, a_left_index) -
          columnCenters(a_spatial_index, 0, a_right_index));
}

void ELVIRA_3D::makeColumnSumsX(void) {
  for (int j = -1; j < 2; ++j) {
    double sum = {0.0};
    for (int i = -1; i < 2; ++i) {
      sum += neighborhood_VF_m->getStoredMoments(i, j, 0) *
             (neighborhood_VF_m->getCell(i, j, 0)).calculateSideLength(0);
    }
    columnCenters(j, 0) =
        (neighborhood_VF_m->getCell(0, j, 0)).calculateCentroid();
    this->columnSum(j, 0) = sum;
  }
  for (int k = -1; k < 2; ++k) {
    double sum = {0.0};
    for (int i = -1; i < 2; ++i) {
      sum += neighborhood_VF_m->getStoredMoments(i, 0, k) *
             (neighborhood_VF_m->getCell(i, 0, k)).calculateSideLength(0);
    }
    columnCenters(0, k) =
        (neighborhood_VF_m->getCell(0, 0, k)).calculateCentroid();
    this->columnSum(0, k) = sum;
  }
}

void ELVIRA_3D::makeColumnSumsY(void) {
  for (int i = -1; i < 2; ++i) {
    double sum = {0.0};
    for (int j = -1; j < 2; ++j) {
      sum += neighborhood_VF_m->getStoredMoments(i, j, 0) *
             (neighborhood_VF_m->getCell(i, j, 0)).calculateSideLength(1);
    }
    columnCenters(i, 0) =
        (neighborhood_VF_m->getCell(i, 0, 0)).calculateCentroid();
    this->columnSum(i, 0) = sum;
  }

  for (int k = -1; k < 2; ++k) {
    double sum = {0.0};
    for (int j = -1; j < 2; ++j) {
      sum += neighborhood_VF_m->getStoredMoments(0, j, k) *
             (neighborhood_VF_m->getCell(0, j, k)).calculateSideLength(1);
    }
    columnCenters(0, k) =
        (neighborhood_VF_m->getCell(0, 0, k)).calculateCentroid();
    this->columnSum(0, k) = sum;
  }
}

void ELVIRA_3D::makeColumnSumsZ(void) {
  for (int i = -1; i < 2; ++i) {
    double sum = {0.0};
    for (int k = -1; k < 2; ++k) {
      sum += neighborhood_VF_m->getStoredMoments(i, 0, k) *
             (neighborhood_VF_m->getCell(i, 0, k)).calculateSideLength(2);
    }
    columnCenters(i, 0) =
        (neighborhood_VF_m->getCell(i, 0, 0)).calculateCentroid();
    this->columnSum(i, 0) = sum;
  }

  for (int j = -1; j < 2; ++j) {
    double sum = {0.0};
    for (int k = -1; k < 2; ++k) {
      sum += neighborhood_VF_m->getStoredMoments(0, j, k) *
             (neighborhood_VF_m->getCell(0, j, k)).calculateSideLength(2);
    }
    columnCenters(0, j) =
        (neighborhood_VF_m->getCell(0, j, 0)).calculateCentroid();
    this->columnSum(0, j) = sum;
  }
}

void ELVIRA_3D::tryNormal(Normal a_normal) {
  a_normal.normalize();
  guess_reconstruction_m[0] = Plane(
      a_normal, findDistanceOnePlane(
                    neighborhood_VF_m->getCell(0, 0, 0),
                    neighborhood_VF_m->getStoredMoments(0, 0, 0), a_normal));
  double error = {0.0};
  for (int k = -1; k < 2; ++k) {
    for (int j = -1; j < 2; ++j) {
      for (int i = -1; i < 2; ++i) {
        double volume_fraction =
            getVolumeFraction<ReconstructionDefaultCuttingMethod>(
                neighborhood_VF_m->getCell(i, j, k), guess_reconstruction_m);
        error += std::pow(
            (volume_fraction - neighborhood_VF_m->getStoredMoments(i, j, k)),
            2);
      }
    }
  }
  if (error < minimum_error_m) {
    best_reconstruction_m = guess_reconstruction_m;
    minimum_error_m = error;
  }
}

//******************************************************************* //
//     Debug function implementations below this.
//******************************************************************* //

PlanarSeparator ELVIRADebug<ELVIRA_2D>::solve(
    const ELVIRANeighborhood* a_neighborhood_pointer) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_VF_m = a_neighborhood_pointer;
  return this->solve();
}

PlanarSeparator ELVIRADebug<ELVIRA_2D>::solve(void) {
  this->writeOutVolumeFractions();
  minimum_error_m = DBL_MAX;
  double tmp_normal_approximations[3];

  // Assume x integration direction
  this->makeColumnSumsX();
  double assumed_normal = neighborhood_VF_m->getStoredMoments(1, 0) >
                                  neighborhood_VF_m->getStoredMoments(-1, 0)
                              ? -1.0
                              : 1.0;
  this->fillOtherNormals(1, tmp_normal_approximations);
  // Try all normal combinations
  for (int y_normal = 0; y_normal < 3; ++y_normal) {
    this->tryNormal(
        Normal(assumed_normal, tmp_normal_approximations[y_normal], 0.0));
  }

  // Assume y integration direction
  this->makeColumnSumsY();
  assumed_normal = neighborhood_VF_m->getStoredMoments(0, 1) >
                           neighborhood_VF_m->getStoredMoments(0, -1)
                       ? -1.0
                       : 1.0;
  this->fillOtherNormals(0, tmp_normal_approximations);

  // Try all normal combinations
  for (int x_normal = 0; x_normal < 3; ++x_normal) {
    this->tryNormal(
        Normal(tmp_normal_approximations[x_normal], assumed_normal, 0.0));
  }

  return best_reconstruction_m;
}

void ELVIRADebug<ELVIRA_2D>::tryNormal(Normal a_normal) {
  this->writeOutPlane(a_normal, "bestPlane", counter_m);
  ++counter_m;
  ELVIRA_2D::tryNormal(a_normal);
}

/// \brief Write out the centroids and weights to
/// to enable visualization of what optimization is driving towards.
void ELVIRADebug<ELVIRA_2D>::writeOutVolumeFractions(void) {
  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      int n = (j + 1) * 3 + (i + 1);
      std::cout << "liq_VF(" << n + 1
                << ") = " << neighborhood_VF_m->getStoredMoments(i, j, -1)
                << "; \n";
    }
  }
}
void ELVIRADebug<ELVIRA_2D>::writeOutPlane(
    Normal a_normal, const std::string& a_prefix,
    const UnsignedIndex_t a_iteration_number) {
  a_normal.normalize();
  PlanarSeparator a_reconstruction = PlanarSeparator::fromOnePlane(Plane(
      a_normal, findDistanceOnePlane(
                    neighborhood_VF_m->getCell(0, 0, -1),
                    neighborhood_VF_m->getStoredMoments(0, 0, -1), a_normal)));
  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      int cell = (j + 1) * 3 + (i + 1);
      for (UnsignedIndex_t plane = 0;
           plane < a_reconstruction.getNumberOfPlanes(); ++plane) {
        Polygon poly = getPlanePolygonFromReconstruction<Polygon>(
            neighborhood_VF_m->getCell(i, j, -1), a_reconstruction,
            a_reconstruction[plane]);
        for (UnsignedIndex_t n = 0; n < poly.getNumberOfVertices(); ++n) {
          std::cout << a_prefix << "(1:3," << n + 1 << "," << cell + 1 << ","
                    << plane + 1 << "," << a_iteration_number + 1 << ") = [";
          std::cout << poly[n].x() << "," << poly[n].y() << "," << poly[n].z()
                    << "];\n";
        }
        std::cout << a_prefix << "Nvert(" << cell + 1 << "," << plane + 1 << ","
                  << a_iteration_number + 1 << ") = ";
        std::cout << poly.getNumberOfVertices() << ";\n";
      }
      std::cout << a_prefix << "NPlane(" << a_iteration_number + 1 << ","
                << cell + 1 << ") = ";
      std::cout << a_reconstruction.getNumberOfPlanes() << ";\n";
    }
  }
}

PlanarSeparator ELVIRADebug<ELVIRA_3D>::solve(
    const ELVIRANeighborhood* a_neighborhood_pointer) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_VF_m = a_neighborhood_pointer;
  return this->solve();
}

PlanarSeparator ELVIRADebug<ELVIRA_3D>::solve(void) {
  this->writeOutVolumeFractions();
  minimum_error_m = DBL_MAX;
  double tmp_normal_approximations[3][2];

  // Assume x integration direction
  this->makeColumnSumsX();
  double assumed_normal = neighborhood_VF_m->getStoredMoments(1, 0, 0) >
                                  neighborhood_VF_m->getStoredMoments(-1, 0, 0)
                              ? -1.0
                              : 1.0;
  this->fillOtherNormals(1, 2, tmp_normal_approximations);
  // Try all normal combinations
  for (int z_normal = 0; z_normal < 3; ++z_normal) {
    for (int y_normal = 0; y_normal < 3; ++y_normal) {
      this->tryNormal(Normal(assumed_normal,
                             tmp_normal_approximations[y_normal][0],
                             tmp_normal_approximations[z_normal][1]));
    }
  }

  // Assume y integration direction
  this->makeColumnSumsY();
  assumed_normal = neighborhood_VF_m->getStoredMoments(0, 1, 0) >
                           neighborhood_VF_m->getStoredMoments(0, -1, 0)
                       ? -1.0
                       : 1.0;
  this->fillOtherNormals(0, 2, tmp_normal_approximations);

  // Try all normal combinations
  for (int z_normal = 0; z_normal < 3; ++z_normal) {
    for (int x_normal = 0; x_normal < 3; ++x_normal) {
      this->tryNormal(Normal(tmp_normal_approximations[x_normal][0],
                             assumed_normal,
                             tmp_normal_approximations[z_normal][1]));
    }
  }

  // Assume z integration direction
  this->makeColumnSumsZ();
  assumed_normal = neighborhood_VF_m->getStoredMoments(0, 0, 1) >
                           neighborhood_VF_m->getStoredMoments(0, 0, -1)
                       ? -1.0
                       : 1.0;
  this->fillOtherNormals(0, 1, tmp_normal_approximations);

  // Try all normal combinations
  for (int y_normal = 0; y_normal < 3; ++y_normal) {
    for (int x_normal = 0; x_normal < 3; ++x_normal) {
      this->tryNormal(Normal(tmp_normal_approximations[x_normal][0],
                             tmp_normal_approximations[y_normal][1],
                             assumed_normal));
    }
  }
  return best_reconstruction_m;
}

void ELVIRADebug<ELVIRA_3D>::tryNormal(Normal a_normal) {
  this->writeOutPlane(a_normal, "bestPlane", counter_m);
  ++counter_m;
  ELVIRA_3D::tryNormal(a_normal);
}

/// \brief Write out the centroids and weights to
/// to enable visualization of what optimization is driving towards.
void ELVIRADebug<ELVIRA_3D>::writeOutVolumeFractions(void) {
  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      int n = (j + 1) * 3 + (i + 1);
      std::cout << "liq_VF(" << n + 1 << ") = "
                << ELVIRA_3D::neighborhood_VF_m->getStoredMoments(i, j, -1)
                << "; \n";
    }
  }
}
void ELVIRADebug<ELVIRA_3D>::writeOutPlane(
    Normal a_normal, const std::string& a_prefix,
    const UnsignedIndex_t a_iteration_number) {
  a_normal.normalize();
  PlanarSeparator a_reconstruction = PlanarSeparator::fromOnePlane(Plane(
      a_normal, findDistanceOnePlane(
                    neighborhood_VF_m->getCell(0, 0, 0),
                    neighborhood_VF_m->getStoredMoments(0, 0, 0), a_normal)));
  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      int cell = (j + 1) * 3 + (i + 1);
      for (UnsignedIndex_t plane = 0;
           plane < a_reconstruction.getNumberOfPlanes(); ++plane) {
        Polygon poly = getPlanePolygonFromReconstruction<Polygon>(
            neighborhood_VF_m->getCell(i, j, -1), a_reconstruction,
            a_reconstruction[plane]);
        for (UnsignedIndex_t n = 0; n < poly.getNumberOfVertices(); ++n) {
          std::cout << a_prefix << "(1:3," << n + 1 << "," << cell + 1 << ","
                    << plane + 1 << "," << a_iteration_number + 1 << ") = [";
          std::cout << poly[n].x() << "," << poly[n].y() << "," << poly[n].z()
                    << "];\n";
        }
        std::cout << a_prefix << "Nvert(" << cell + 1 << "," << plane + 1 << ","
                  << a_iteration_number + 1 << ") = ";
        std::cout << poly.getNumberOfVertices() << ";\n";
      }
      std::cout << a_prefix << "NPlane(" << a_iteration_number + 1 << ","
                << cell + 1 << ") = ";
      std::cout << a_reconstruction.getNumberOfPlanes() << ";\n";
    }
  }
}

}  // namespace IRL
