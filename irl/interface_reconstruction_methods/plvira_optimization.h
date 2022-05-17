// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_OPTIMIZATION_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_OPTIMIZATION_H_

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>  // Eigen header

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/interface_reconstruction_methods/optimization_behavior.h"
#include "irl/interface_reconstruction_methods/plvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/reconstruction_cleaning.h"
#include "irl/optimization/bfgs.h"
#include "irl/optimization/constrained_levenberg_marquardt.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/compiler_type.h"
#include "irl/parameters/defined_types.h"
#include "irl/planar_reconstruction/planar_separator.h"

namespace IRL {

// static constexpr UnsignedIndex_t PLVIRA_2D_columns = 1;
static constexpr UnsignedIndex_t PLVIRA_3D_parameters = 6;
static constexpr UnsignedIndex_t PLVIRA_3D_constraints = 1;

// template <class CellType>
// class PLVIRA_2D;
template <class CellType>
class PLVIRA_3D;

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
class PLVIRACommon {
  //   friend PLVIRA_2D<CellType>;
  friend PLVIRA_3D<CellType>;

 public:
  PLVIRACommon(void) = default;

  template <class PLVIRAType>
  Paraboloid runOptimization(
      PLVIRAType* a_ptr_to_PLVIRA_object,
      const PLVIRANeighborhood<CellType>& a_neighborhood_geometry,
      const Paraboloid& a_reconstruction);

  Paraboloid getFinalReconstruction(void);

  void setOptimizationBehavior(const OptimizationBehavior& a_parameters);

  Eigen::Matrix<double, Eigen::Dynamic, 1> calculateVectorError(void);

  bool errorTooHigh(const double a_error);

  bool iterationTooHigh(const UnsignedIndex_t a_iteration);

  void increaseLambda(double* a_lambda);

  void decreaseLambda(double* a_lambda);

  bool shouldComputeJacobian(const UnsignedIndex_t a_iteration,
                             const UnsignedIndex_t a_last_jacobian);

#ifndef USING_INTEL_COMPILER
  bool minimumReached(
      const Eigen::Matrix<double, static_cast<int>(kParameters + kConstraints),
                          1>& a_delta) const;

  Eigen::Matrix<double, static_cast<int>(kParameters + kConstraints), 1>
  getJacobianStepSize(void) const;
#else  // Intel really doesn't like the static cast of kColumns
  bool minimumReached(const Eigen::Matrix<double, kParameters + kConstraints,
                                          1>& a_delta) const;

  Eigen::Matrix<double, kParameters + kConstraints, 1> getJacobianStepSize(
      void) const;
#endif

  double calculateScalarError(void);

  void updateBestGuess(void);

  Eigen::Matrix<double, Eigen::Dynamic, 1> calculateChangeInGuess(void);

  const Paraboloid& getBestReconstruction(void);

  const Paraboloid& getGuessReconstruction(void);

  const ReferenceFrame& getBestReferenceFrame(void);

  const ReferenceFrame& getGuessReferenceFrame(void);

  void allocateMatrices(const UnsignedIndex_t a_neighborhood_size);

  void fillGeometryAndWeightVectors(void);

  void setWeightedGeometryVectorFromReconstruction(
      const Paraboloid& a_reconstruction);

  //   Paraboloid getReconstructionFromLVIRAParam(
  //       const Paraboloid& a_reference_frame);

  ~PLVIRACommon(void) = default;

 public:
  const PLVIRANeighborhood<CellType>* neighborhood_m;
  Eigen::Matrix<double, Eigen::Dynamic, 1> weights_m;
  Eigen::Matrix<double, Eigen::Dynamic, 1> correct_values_m;
  OptimizationBehavior optimization_behavior_m;
  Eigen::Matrix<double, Eigen::Dynamic, 1> guess_values_m;
  Paraboloid guess_reconstruction_m;
  ReferenceFrame guess_reference_frame_m;
  Eigen::Matrix<double, Eigen::Dynamic, 1> best_values_m;
  Paraboloid best_reconstruction_m;
  ReferenceFrame best_reference_frame_m;
};

template <class CellType>
class PLVIRA_3D : public PLVIRACommon<CellType, PLVIRA_3D_parameters,
                                      PLVIRA_3D_constraints> {
 public:
  using cell_type = CellType;
  static constexpr UnsignedIndex_t parameters_m = PLVIRA_3D_parameters;
  static constexpr UnsignedIndex_t constraints_m = PLVIRA_3D_constraints;
  static constexpr UnsignedIndex_t columns_m =
      PLVIRA_3D_parameters + PLVIRA_3D_constraints;

  PLVIRA_3D(void) = default;
  ~PLVIRA_3D(void) = default;

  Paraboloid solve(const PLVIRANeighborhood<CellType>& a_neighborhood,
                   const Paraboloid& a_reconstruction);

  void setup(const Paraboloid& a_reconstruction);

  void updateGuess(const Eigen::Matrix<double, columns_m, 1>* const a_delta);

 private:
  UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);
};

}  // namespace IRL

#include "irl/interface_reconstruction_methods/plvira_optimization.tpp"

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_OPTIMIZATION_H_
