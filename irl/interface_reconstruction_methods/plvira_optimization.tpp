// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_OPTIMIZATION_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_OPTIMIZATION_TPP_

namespace IRL {

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
template <class PLVIRAType>
Paraboloid PLVIRACommon<CellType, kParameters, kConstraints>::runOptimization(
    PLVIRAType* a_ptr_to_PLVIRA_object,
    const PLVIRANeighborhood<CellType>& a_neighborhood_geometry,
    const Paraboloid& a_reconstruction) {
  neighborhood_m = &a_neighborhood_geometry;
  a_ptr_to_PLVIRA_object->setup(a_reconstruction);
  ConstrainedLevenbergMarquardt<PLVIRAType, -1,
                                static_cast<int>(PLVIRAType::parameters_m),
                                static_cast<int>(PLVIRAType::constraints_m)>
      lm_solver;
  lm_solver.solve(a_ptr_to_PLVIRA_object);
  return a_ptr_to_PLVIRA_object->getFinalReconstruction();
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
Paraboloid PLVIRACommon<CellType, kParameters,
                        kConstraints>::getFinalReconstruction(void) {
  return this->getBestReconstruction();
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
void PLVIRACommon<CellType, kParameters, kConstraints>::setOptimizationBehavior(
    const OptimizationBehavior& a_parameters) {
  optimization_behavior_m = a_parameters;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
Eigen::Matrix<double, Eigen::Dynamic, 1>
PLVIRACommon<CellType, kParameters, kConstraints>::calculateVectorError(void) {
  return (correct_values_m - guess_values_m);
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
bool PLVIRACommon<CellType, kParameters, kConstraints>::errorTooHigh(
    const double a_error) {
  return a_error > optimization_behavior_m.acceptable_error;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
bool PLVIRACommon<CellType, kParameters, kConstraints>::iterationTooHigh(
    const UnsignedIndex_t a_iteration) {
  return a_iteration > optimization_behavior_m.maximum_iterations;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
void PLVIRACommon<CellType, kParameters, kConstraints>::increaseLambda(
    double* a_lambda) {
  (*a_lambda) *= optimization_behavior_m.lambda_increase;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
void PLVIRACommon<CellType, kColumns>::decreaseLambda(double* a_lambda) {
  (*a_lambda) *= optimization_behavior_m.lambda_decrease;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
bool PLVIRACommon<CellType, kParameters, kConstraints>::shouldComputeJacobian(
    const UnsignedIndex_t a_iteration, const UnsignedIndex_t a_last_jacobian) {
  return a_iteration - a_last_jacobian >
         optimization_behavior_m.delay_jacobian_amount;
}

#ifndef USING_INTEL_COMPILER
template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
bool PLVIRACommon<CellType, kParameters, kConstraints>::minimumReached(
    const Eigen::Matrix<double, static_cast<int>(kParameters + kConstraints),
                        1>& a_delta) const {
  return a_delta.squaredNorm() <
         PLVIRACommon<CellType, kParameters,
                      kConstraints>::optimization_behavior_m
                 .minimum_angle_change *
             PLVIRACommon<CellType, kParameters,
                          kConstraints>::optimization_behavior_m
                 .minimum_angle_change;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
Eigen::Matrix<double, static_cast<int>(kParameters + kConstraints), 1>
PLVIRACommon<CellType, kParameters, kConstraints>::getJacobianStepSize(
    void) const {
  const double angle_to_use =
      PLVIRACommon<CellType, kParameters, kConstraints>::optimization_behavior_m
          .finite_difference_angle;
  return Eigen::Matrix<double, static_cast<int>(kParameters + kConstraints),
                       1>::Constant(angle_to_use);
}

#else
template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
bool PLVIRACommon<CellType, kParameters, kConstraints>::minimumReached(
    const Eigen::Matrix<double, kParameters + kConstraints, 1>& a_delta) const {
  return a_delta.squaredNorm() <
         PLVIRACommon<CellType, kParameters,
                      kConstraints>::optimization_behavior_m
                 .minimum_angle_change *
             PLVIRACommon<CellType, kParameters,
                          kConstraints>::optimization_behavior_m
                 .minimum_angle_change;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
Eigen::Matrix<double, kParameters + kConstraints, 1>
PLVIRACommon<CellType, kParameters, kConstraints>::getJacobianStepSize(
    void) const {
  const double angle_to_use =
      PLVIRACommon<CellType, kParameters, kConstraints>::optimization_behavior_m
          .finite_difference_angle;
  return Eigen::Matrix<double, static_cast<int>(kParameters + kConstraints),
                       1>::Constant(angle_to_use);
}
#endif

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
double PLVIRACommon<CellType, kParameters, kConstraints>::calculateScalarError(
    void) {
  return (correct_values_m - guess_values_m).squaredNorm();
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
void PLVIRACommon<CellType, kParameters, kConstraints>::updateBestGuess(void) {
  best_values_m = guess_values_m;
  best_reconstruction_m = guess_reconstruction_m;
  best_reference_frame_m = guess_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
Eigen::Matrix<double, Eigen::Dynamic, 1>
PLVIRACommon<CellType, kParameters, kConstraints>::calculateChangeInGuess(
    void) {
  return guess_values_m - best_values_m;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
const Paraboloid&
PLVIRACommon<CellType, kParameters, kConstraints>::getBestReconstruction(void) {
  return best_reconstruction_m;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
const Paraboloid& PLVIRACommon<CellType, kParameters,
                               kConstraints>::getGuessReconstruction(void) {
  return guess_reconstruction_m;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
const ReferenceFrame&
PLVIRACommon<CellType, kParameters, kConstraints>::getBestReferenceFrame(void) {
  return best_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
const ReferenceFrame& PLVIRACommon<CellType, kParameters,
                                   kConstraints>::getGuessReferenceFrame(void) {
  return guess_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
void PLVIRACommon<CellType, kParameters, kConstraints>::allocateMatrices(
    const UnsignedIndex_t a_neighborhood_size) {
  const int rows = static_cast<int>(a_neighborhood_size);
  this->weights_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->correct_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->guess_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->best_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
}

template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
void PLVIRACommon<CellType, kParameters, kConstraints>::
    setWeightedGeometryVectorFromReconstruction(
        const Paraboloid& a_reconstruction) {
  for (UnsignedIndex_t i = 0; i < neighborhood_m->size(); ++i) {
    const double volume_fraction =
        getVolumeFraction<ReconstructionDefaultCuttingMethod>(
            neighborhood_m->getCell(i), a_reconstruction);
    guess_values_m(i) = weights_m(i) * volume_fraction;
  }
}

// Turn off warnings about sign conversion because need to work
// with Eigen which using long int
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
template <class CellType, UnsignedIndex_t kParameters,
          UnsignedIndex_t kConstraints>
void PLVIRACommon<CellType, kParameters,
                  kConstraints>::fillGeometryAndWeightVectors(void) {
  for (UnsignedIndex_t n = 0; n < neighborhood_m->size(); ++n) {
    // Add correct volume fractions
    correct_values_m(n) = neighborhood_m->getStoredMoments(n);

    // Add even weighting of 1 for liquid volume
    weights_m(n) = 1.0;
  }

  // Normalize weight to have an L2 magnitude of 1
  const double liquid_volume_weight_magnitude = 1.0 / weights_m.norm();
  for (UnsignedIndex_t n = 0; n < guess_values_m.rows(); ++n) {
    weights_m(n) *= liquid_volume_weight_magnitude;
    correct_values_m(n) *= weights_m(n);
  }
}
#pragma GCC diagnostic pop

// template <class CellType, UnsignedIndex_t kColumns>
// Paraboloid
// PLVIRACommon<CellType, kColumns>::getReconstructionFromLVIRAParam(
//     const ReferenceFrame& a_reference_frame) {
//   auto planar_separator = PlanarSeparator::fromOnePlane(
//       Plane(a_reference_frame[2],
//             a_reference_frame[2] *
//                 this->neighborhood_m->getCenterCell().calculateCentroid()));
//   setDistanceToMatchVolumeFractionPartialFill(
//       this->neighborhood_m->getCenterCell(),
//       this->neighborhood_m->getCenterCellStoredMoments(), &planar_separator);
//   return planar_separator;
// }

//******************************************************************* //
//     Functions for PLVIRA_3D below this
//******************************************************************* //
template <class CellType>
Paraboloid PLVIRA_3D<CellType>::solve(
    const PLVIRANeighborhood<CellType>& a_neighborhood,
    const Paraboloid& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class CellType>
void PLVIRA_3D<CellType>::setup(const Paraboloid& a_reconstruction) {
  //   this->allocateMatrices(this->neighborhood_m->size());
  //   this->fillGeometryAndWeightVectors();
  //   this->best_reference_frame_m =
  //       getOrthonormalSystem(a_reconstruction[0].normal());
}

template <class CellType>
void PLVIRA_3D<CellType>::updateGuess(
    const Eigen::Matrix<double, columns_m, 1>* const a_delta) {
  //   UnitQuaternion rotation =
  //       this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  //   this->guess_reference_frame_m = rotation * this->best_reference_frame_m;
  //   this->guess_reconstruction_m =
  //       this->getReconstructionFromLVIRAParam(this->guess_reference_frame_m);
  //   this->setWeightedGeometryVectorFromReconstruction(
  //       this->guess_reconstruction_m);
}

template <class CellType>
UnitQuaternion PLVIRA_3D<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(1), a_reference_frame[1]) *
         UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_OPTIMIZATION_TPP_
