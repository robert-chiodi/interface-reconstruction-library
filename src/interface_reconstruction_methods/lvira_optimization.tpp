// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_OPTIMIZATION_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_OPTIMIZATION_TPP_

namespace IRL {

template <class CellType, UnsignedIndex_t kColumns>
template <class LVIRAType>
PlanarSeparator LVIRACommon<CellType, kColumns>::runOptimization(
    LVIRAType* a_ptr_to_LVIRA_object,
    const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
    const PlanarSeparator& a_reconstruction) {
  neighborhood_m = &a_neighborhood_geometry;
  a_ptr_to_LVIRA_object->setup(a_reconstruction);
  LevenbergMarquardt<LVIRAType, -1, static_cast<int>(LVIRAType::columns_m)>
      lm_solver;
  lm_solver.solve(a_ptr_to_LVIRA_object,
                  static_cast<int>(correct_values_m.rows()),
                  a_ptr_to_LVIRA_object->getJacobianStepSize());
  return a_ptr_to_LVIRA_object->getFinalReconstruction();
}

template <class CellType, UnsignedIndex_t kColumns>
PlanarSeparator LVIRACommon<CellType, kColumns>::getFinalReconstruction(void) {
  return this->getBestReconstruction();
}

template <class CellType, UnsignedIndex_t kColumns>
Eigen::Matrix<double, Eigen::Dynamic, 1>
LVIRACommon<CellType, kColumns>::calculateVectorError(void) {
  return (correct_values_m - guess_values_m);
}

template <class CellType, UnsignedIndex_t kColumns>
bool LVIRACommon<CellType, kColumns>::errorTooHigh(const double a_error) {
  return a_error > acceptable_error_m;
}

template <class CellType, UnsignedIndex_t kColumns>
bool LVIRACommon<CellType, kColumns>::iterationTooHigh(
    const UnsignedIndex_t a_iteration) {
  return a_iteration > maximum_iterations_m;
}

template <class CellType, UnsignedIndex_t kColumns>
void LVIRACommon<CellType, kColumns>::increaseLambda(double* a_lambda) {
  (*a_lambda) *= lambda_increase_m;
}

template <class CellType, UnsignedIndex_t kColumns>
void LVIRACommon<CellType, kColumns>::decreaseLambda(double* a_lambda) {
  (*a_lambda) *= lambda_decrease_m;
}

template <class CellType, UnsignedIndex_t kColumns>
bool LVIRACommon<CellType, kColumns>::shouldComputeJacobian(
    const UnsignedIndex_t a_iteration, const UnsignedIndex_t a_last_jacobian) {
  return a_iteration - a_last_jacobian > delay_jacobian_amount_m;
}

#ifndef USING_INTEL_COMPILER
template <class CellType, UnsignedIndex_t kColumns>
bool LVIRACommon<CellType, kColumns>::minimumReached(
    const Eigen::Matrix<double, static_cast<int>(kColumns), 1>& a_delta) const {
  return a_delta.squaredNorm() <
         LVIRACommon<CellType, kColumns>::minimum_angle_change_m *
             LVIRACommon<CellType, kColumns>::minimum_angle_change_m;
}

template <class CellType, UnsignedIndex_t kColumns>
Eigen::Matrix<double, static_cast<int>(kColumns), 1>
LVIRACommon<CellType, kColumns>::getJacobianStepSize(void) const {
  const double angle_to_use =
      LVIRACommon<CellType, kColumns>::fininite_difference_angle_m;
  return Eigen::Matrix<double, static_cast<int>(kColumns), 1>::Constant(
      angle_to_use);
}

#else
template <class CellType, UnsignedIndex_t kColumns>
bool LVIRACommon<CellType, kColumns>::minimumReached(
    const Eigen::Matrix<double, kColumns, 1>& a_delta) const {
  return a_delta.squaredNorm() <
         LVIRACommon<CellType, kColumns>::minimum_angle_change_m *
             LVIRACommon<CellType, kColumns>::minimum_angle_change_m;
}

template <class CellType, UnsignedIndex_t kColumns>
Eigen::Matrix<double, kColumns, 1>
LVIRACommon<CellType, kColumns>::getJacobianStepSize(void) const {
  const double angle_to_use =
      LVIRACommon<CellType, kColumns>::fininite_difference_angle_m;
  return Eigen::Matrix<double, static_cast<int>(kColumns), 1>::Constant(
      angle_to_use);
}
#endif

template <class CellType, UnsignedIndex_t kColumns>
double LVIRACommon<CellType, kColumns>::calculateScalarError(void) {
  return (correct_values_m - guess_values_m).squaredNorm();
}

template <class CellType, UnsignedIndex_t kColumns>
void LVIRACommon<CellType, kColumns>::updateBestGuess(void) {
  best_values_m = guess_values_m;
  best_reconstruction_m = guess_reconstruction_m;
  best_reference_frame_m = guess_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kColumns>
Eigen::Matrix<double, Eigen::Dynamic, 1>
LVIRACommon<CellType, kColumns>::calculateChangeInGuess(void) {
  return guess_values_m - best_values_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const PlanarSeparator& LVIRACommon<CellType, kColumns>::getBestReconstruction(
    void) {
  return best_reconstruction_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const PlanarSeparator& LVIRACommon<CellType, kColumns>::getGuessReconstruction(
    void) {
  return guess_reconstruction_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const ReferenceFrame& LVIRACommon<CellType, kColumns>::getBestReferenceFrame(
    void) {
  return best_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const ReferenceFrame& LVIRACommon<CellType, kColumns>::getGuessReferenceFrame(
    void) {
  return guess_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kColumns>
void LVIRACommon<CellType, kColumns>::allocateMatrices(
    const UnsignedIndex_t a_neighborhood_size) {
  const int rows = static_cast<int>(a_neighborhood_size);
  this->weights_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->correct_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->guess_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->best_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
}

template <class CellType, UnsignedIndex_t kColumns>
void LVIRACommon<CellType, kColumns>::
    setWeightedGeometryVectorFromReconstruction(
        const PlanarSeparator& a_reconstruction) {
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
template <class CellType, UnsignedIndex_t kColumns>
void LVIRACommon<CellType, kColumns>::fillGeometryAndWeightVectors(void) {
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

template <class CellType, UnsignedIndex_t kColumns>
PlanarSeparator
LVIRACommon<CellType, kColumns>::getReconstructionFromLVIRAParam(
    const ReferenceFrame& a_reference_frame) {
  auto planar_separator = PlanarSeparator::fromOnePlane(
      Plane(a_reference_frame[2],
            a_reference_frame[2] *
                this->neighborhood_m->getCenterCell().calculateCentroid()));
  setDistanceToMatchVolumeFractionPartialFill(
      this->neighborhood_m->getCenterCell(),
      this->neighborhood_m->getCenterCellStoredMoments(), &planar_separator);
  return planar_separator;
}

//******************************************************************* //
//     Functions for LVIRA_2D below this
//******************************************************************* //
/// \brief Perform optimization and return planar separator.
template <class CellType>
PlanarSeparator LVIRA_2D<CellType>::solve(
    const LVIRANeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class CellType>
void LVIRA_2D<CellType>::setup(const PlanarSeparator& a_reconstruction) {
  this->allocateMatrices(this->neighborhood_m->size());
  this->fillGeometryAndWeightVectors();
  Normal tangent0 = Normal(0.0, 0.0, 1.0);
  Normal tangent1 =
      crossProductNormalized(tangent0, a_reconstruction[0].normal());
  this->best_reference_frame_m =
      ReferenceFrame(tangent0, tangent1, a_reconstruction[0].normal());
}

template <class CellType>
void LVIRA_2D<CellType>::updateGuess(
    const Eigen::Matrix<double, columns_m, 1>* const a_delta) {
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * this->best_reference_frame_m;
  this->guess_reconstruction_m =
      this->getReconstructionFromLVIRAParam(this->guess_reference_frame_m);
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
}

template <class CellType>
UnitQuaternion LVIRA_2D<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

//******************************************************************* //
//     Functions for LVIRA_3D below this
//******************************************************************* //
template <class CellType>
PlanarSeparator LVIRA_3D<CellType>::solve(
    const LVIRANeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class CellType>
void LVIRA_3D<CellType>::setup(const PlanarSeparator& a_reconstruction) {
  this->allocateMatrices(this->neighborhood_m->size());
  this->fillGeometryAndWeightVectors();
  this->best_reference_frame_m =
      getOrthonormalSystem(a_reconstruction[0].normal());
}

template <class CellType>
void LVIRA_3D<CellType>::updateGuess(
    const Eigen::Matrix<double, columns_m, 1>* const a_delta) {
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * this->best_reference_frame_m;
  this->guess_reconstruction_m =
      this->getReconstructionFromLVIRAParam(this->guess_reference_frame_m);
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
}

template <class CellType>
UnitQuaternion LVIRA_3D<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(1), a_reference_frame[1]) *
         UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

//******************************************************************* //
//     Functions for Debug LVIRA versions below this
//******************************************************************* //

template <class LVIRAType>
void LVIRADebug<LVIRAType>::updateGuess(
    const Eigen::Matrix<double, LVIRAType::columns_m, 1>* const a_delta) {
  LVIRAType::updateGuess(a_delta);
  guess_reconstruction_history.push_back(LVIRAType::getGuessReconstruction());
}

template <class LVIRAType>
PlanarSeparator LVIRADebug<LVIRAType>::solve(
    const LVIRANeighborhood<typename LVIRAType::cell_type>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class LVIRAType>
void LVIRADebug<LVIRAType>::updateBestGuess(void) {
  LVIRAType::updateBestGuess();
  best_reconstruction_history.push_back(LVIRAType::getBestReconstruction());
  PlanarSeparator tmp = best_reconstruction_history.back();
  std::cout << "**************************************** \n";
  std::cout << "Writing out plane: " << best_reconstruction_history.size() - 1
            << '\n';
  std::cout << tmp << std::endl;
}

template <class LVIRAType>
PlanarSeparator LVIRADebug<LVIRAType>::getFinalReconstruction(void) {
  this->writeOutVolumesAndWeights();
  for (std::size_t n = 0; n < best_reconstruction_history.size(); ++n) {
    this->writeOutPlane(best_reconstruction_history[n], "bestPlane", n);
  }
  return LVIRAType::getFinalReconstruction();
}

template <class LVIRAType>
void LVIRADebug<LVIRAType>::writeOutPlane(
    const PlanarSeparator& a_reconstruction, const std::string& a_prefix,
    const std::size_t a_iteration_number) {
  for (UnsignedIndex_t cell = 0; cell < this->neighborhood_m->size(); ++cell) {
    for (UnsignedIndex_t plane = 0;
         plane < a_reconstruction.getNumberOfPlanes(); ++plane) {
      Polygon poly = getPlanePolygonFromReconstruction<Polygon>(
          this->neighborhood_m->getCell(cell), a_reconstruction,
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

template <class LVIRAType>
void LVIRADebug<LVIRAType>::writeOutVolumesAndWeights(void) {
  for (UnsignedIndex_t n = 0; n < this->neighborhood_m->size(); ++n) {
    std::cout << "Liquid_VF(1:4," << n + 1 << ") = [";
    std::cout << this->neighborhood_m->getStoredMoments(n) << " "
              << this->weights_m(n) << "];\n";
  }
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_OPTIMIZATION_TPP_
