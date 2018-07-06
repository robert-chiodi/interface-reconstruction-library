// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_ELVIRA_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_ELVIRA_H_

#include <float.h>

#include <cassert>
#include <string>

#include "src/generic_cutting/cut_polygon.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polygons/polygon.h"
#include "src/interface_reconstruction_methods/elvira_neighborhood.h"
#include "src/moments/cell_collection.h"
#include "src/moments/cell_grouped_moments.h"
#include "src/parameters/defined_types.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

// Forward declare to friend with ELVIRA_2D and ELVIRA_3D.
template <class BaseELVIRA>
class ELVIRADebug;

class ELVIRA_2D {
  friend ELVIRADebug<ELVIRA_2D>;

 public:
  /// \brief Default constructor.
  ELVIRA_2D(void) = default;

  /// \brief Solve the system for the reconstruction, restarting
  /// the neighboring geoemtry
  PlanarSeparator solve(const ELVIRANeighborhood* a_neighborhood_pointer);

  /// \brief Default destructor.
  ~ELVIRA_2D(void) = default;

 private:
  /// \brief Solve the system for the reconstruction.
  PlanarSeparator solve(void);

  /// \brief Return value of column sum at index.
  /// range for index of -1:1, with center cell column
  /// being at index 0,0.
  double& columnSum(const int i);

  /// \brief Return point of column center at index i,k,
  /// with the integrated dimensions index being 0.
  /// range for index of -1:1, with center cell column
  /// being at index 0,0.
  Pt& columnCenters(const int i);

  /// \brief Overloaded column centers to return the location
  /// of the column center for the `dimension`.
  double columnCenters(const UnsignedIndex_t dimension, const int i);

  UnsignedIndex_t calculateLinearIndex(const int i);

  void fillOtherNormals(const UnsignedIndex_t a_dimension_index,
                        double a_normal_approximation_holder[3]);

  /// \brief Make column sums from summing in the x direction.
  void makeColumnSumsX(void);

  /// \brief Make column sums from summing in the y direction.
  void makeColumnSumsY(void);

  /// \brief Compute derivative of the column heights
  double computeDerivative0(const int a_left_index, const int a_right_index,
                            const UnsignedIndex_t spatial_index);

  /// \brief Calculate error given the supplied normal.
  /// If lowest error yet, save as best_reconstruction_m.
  void tryNormal(Normal a_normal);
  /// \brief Storage of the stencil information
  const ELVIRANeighborhood* neighborhood_VF_m;
  /// \brief Array of column sums needed in ELVIRA
  std::array<double, 3> column_sums_m;
  /// \brief Array of points for centers of columns.
  /// Needed for derivatives.
  std::array<Pt, 3> column_centers_m;
  /// \brief Best reconstruction found so far.
  PlanarSeparator guess_reconstruction_m;
  /// \brief Best reconstruction found so far.
  PlanarSeparator best_reconstruction_m;
  /// \brief Lowest error resulting from best reconstruction.
  double minimum_error_m;
};

class ELVIRA_3D {
  friend ELVIRADebug<ELVIRA_3D>;

 public:
  /// \brief Default constructor.
  ELVIRA_3D(void) = default;

  /// \brief Solve the system for the reconstruction, restarting
  /// the neighboring geoemtry
  PlanarSeparator solve(const ELVIRANeighborhood* a_neighborhood_pointer);

  /// \brief Default destructor.
  ~ELVIRA_3D(void) = default;

 private:
  /// \brief Solve the system for the reconstruction.
  PlanarSeparator solve(void);

  /// \brief Return value of column sum at index.
  /// range for index of -1:1, with center cell column
  /// being at index 0,0.
  double& columnSum(const int i, const int j);

  /// \brief Return point of column center at index i,k,
  /// with the integrated dimensions index being 0.
  /// range for index of -1:1, with center cell column
  /// being at index 0,0.
  Pt& columnCenters(const int i, const int j);

  /// \brief Overloaded column centers to return the location
  /// of the column center for the `dimension`.
  double columnCenters(const UnsignedIndex_t dimension, const int i,
                       const int j);

  UnsignedIndex_t calculateLinearIndex(const int i, const int j);

  /// \brief Fill in normal approximations from derivatives of column heights.
  void fillOtherNormals(const UnsignedIndex_t a_first_column_index,
                        const UnsignedIndex_t a_second_column_index,
                        double a_normal_approximation_holder[3][2]);

  /// \brief Make column sums from summing in the x direction.
  void makeColumnSumsX(void);

  /// \brief Make column sums from summing in the y direction.
  void makeColumnSumsY(void);

  /// \brief Make column sums from summing in the z direction.
  void makeColumnSumsZ(void);

  /// \brief Compute derivative of first dimension of column heights.
  double computeDerivative0(const int a_left_index, const int a_right_index,
                            const UnsignedIndex_t a_spatial_index);

  /// \brief Compute derivative of second dimension of column heights.
  double computeDerivative1(const int a_left_index, const int a_right_index,
                            const UnsignedIndex_t a_spatial_index);

  /// \brief Calculate error given the supplied normal.
  /// If lowest error yet, save as best_reconstruction_m.
  void tryNormal(Normal a_normal);

  /// \brief Storage of the stencil information
  const ELVIRANeighborhood* neighborhood_VF_m;
  /// \brief Array of column sums needed in ELVIRA
  std::array<double, 9> column_sums_m;
  /// \brief Array of points for centers of columns.
  /// Needed for derivatives.
  std::array<Pt, 9> column_centers_m;
  /// \brief Best reconstruction found so far.
  PlanarSeparator guess_reconstruction_m;
  /// \brief Best reconstruction found so far.
  PlanarSeparator best_reconstruction_m;
  /// \brief Lowest error resulting from best reconstruction.
  double minimum_error_m;
};

template <>
class ELVIRADebug<ELVIRA_2D> : private ELVIRA_2D {
 public:
  /// \brief Default constructor.
  ELVIRADebug(void) = default;

  /// \brief Solve the system for the reconstruction, restarting
  /// the neighboring geoemtry
  PlanarSeparator solve(const ELVIRANeighborhood* a_neighborhood_pointer);

 private:
  /// \brief Solve the system for the reconstruction.
  PlanarSeparator solve(void);

  void writeOutVolumeFractions(void);

  void writeOutPlane(Normal a_normal, const std::string& a_prefix,
                     const UnsignedIndex_t a_iteration_number);

  void tryNormal(Normal a_normal);

  /// \brief Counter to note which attempted normal it is.
  UnsignedIndex_t counter_m = 0;
};

template <>
class ELVIRADebug<ELVIRA_3D> : private ELVIRA_3D {
 public:
  /// \brief Default constructor.
  ELVIRADebug(void) = default;

  /// \brief Solve the system for the reconstruction, restarting
  /// the neighboring geoemtry
  PlanarSeparator solve(const ELVIRANeighborhood* a_neighborhood_pointer);

 private:
  /// \brief Solve the system for the reconstruction.
  PlanarSeparator solve(void);

  void writeOutVolumeFractions(void);

  void writeOutPlane(Normal a_normal, const std::string& a_prefix,
                     const UnsignedIndex_t a_iteration_number);

  void tryNormal(Normal a_normal);

  /// \brief Counter to note which attempted normal it is.
  UnsignedIndex_t counter_m = 0;
};

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_ELVIRA_H_
