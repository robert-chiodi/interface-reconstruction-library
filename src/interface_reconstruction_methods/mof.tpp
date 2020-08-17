// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_MOF_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_MOF_TPP_

namespace IRL {

template <class CellType>
template <class MOFType>
PlanarSeparator MOFCommon<CellType>::runOptimization(
    MOFType* a_ptr_to_MOF_object,
    const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_moments,
    const double a_internal_weight, const double a_external_weight) {
  assert(a_ptr_to_MOF_object != nullptr);
  a_ptr_to_MOF_object->setup(a_cell_grouped_moments, a_internal_weight,
                             a_external_weight);
  LevenbergMarquardt<MOFType, MOFType::rows_m, MOFType::columns_m> lm_solver;
  lm_solver.solve(a_ptr_to_MOF_object,
                  a_ptr_to_MOF_object->getDefaultInitialDelta());
  //  BFGS<MOFType, MOFType::columns_m> bfgs_solver;
  //  bfgs_solver.solve(a_ptr_to_MOF_object,
  //                    a_ptr_to_MOF_object->getDefaultInitialDelta());
  return a_ptr_to_MOF_object->getFinalReconstruction();
}

template <class CellType>
PlanarSeparator MOFCommon<CellType>::getFinalReconstruction(void) {
  return this->getBestReconstruction();
}

template <class CellType>
double MOFCommon<CellType>::calculateScalarError(void) {
  return (correct_values_m - guess_values_m).squaredNorm();
}

template <class CellType>
Eigen::Matrix<double, 6, 1> MOFCommon<CellType>::calculateVectorError(void) {
  return correct_values_m - guess_values_m;
}

template <class CellType>
Eigen::Matrix<double, 6, 1> MOFCommon<CellType>::calculateChangeInGuess(void) {
  return guess_values_m - best_values_m;
}

template <class CellType>
void MOFCommon<CellType>::increaseLambda(double* a_lambda) {
  (*a_lambda) *= lambda_increase_m;
}

template <class CellType>
void MOFCommon<CellType>::decreaseLambda(double* a_lambda) {
  (*a_lambda) *= lambda_decrease_m;
}

template <class CellType>
bool MOFCommon<CellType>::errorTooHigh(const double a_error) {
  return a_error > acceptable_error_m;
}

template <class CellType>
bool MOFCommon<CellType>::iterationTooHigh(const UnsignedIndex_t a_iteration) {
  return a_iteration > maximum_iterations_m;
}

template <class CellType>
bool MOFCommon<CellType>::shouldComputeJacobian(
    const UnsignedIndex_t a_iteration, const UnsignedIndex_t a_last_jacobian) {
  return a_iteration - a_last_jacobian > delay_jacobian_amount_m;
}

template <class CellType>
void MOFCommon<CellType>::updateBestGuess(void) {
  best_values_m = guess_values_m;
  best_reconstruction_m = guess_reconstruction_m;
  best_reference_frame_m = guess_reference_frame_m;
}

template <class CellType>
void MOFCommon<CellType>::fillGeometryAndWeightVectors(double a_liquid_weight,
                                                       double a_gas_weight) {
  // Set scaling factor to vol^(1/3)
  double encompassing_volume = cell_grouped_data_m->getCell().calculateVolume();
  geom_scale_factor_m = std::pow(encompassing_volume, 1.0 / 3.0);
  volume_fraction_m = this->getInternalVolume() / encompassing_volume;

  // Setup weights and correct_values vector
  double sum = a_liquid_weight + a_gas_weight;
  a_liquid_weight /= (sum * geom_scale_factor_m);
  a_gas_weight /= (sum * geom_scale_factor_m);
  for (UnsignedIndex_t n = 0; n < 3; ++n) {
    weights_m(n) = a_liquid_weight;
    weights_m(3 + n) = a_gas_weight;
    correct_values_m(n) = this->getInternalCentroid()[n];
    correct_values_m(3 + n) = this->getExternalCentroid()[n];
  }
  correct_values_m = correct_values_m.cwiseProduct(weights_m);
}

template <class CellType>
inline PlanarSeparator MOFCommon<CellType>::getReconstructionFromParam(
    const ReferenceFrame& a_reference_frame) {
  const double initial_guess_dist =
      a_reference_frame[2] *
      cell_grouped_data_m->getStoredMoments()[0].centroid();
  PlanarSeparator separator_to_return = PlanarSeparator::fromOnePlane(
      Plane(a_reference_frame[2], initial_guess_dist));
  setDistanceToMatchVolumeFractionPartialFill(
      cell_grouped_data_m->getCell(), volume_fraction_m, &separator_to_return);
  return separator_to_return;
}

template <class CellType>
inline void MOFCommon<CellType>::setWeightedGeometryVectorFromReconstruction(
    const PlanarSeparator& a_reconstruction) {
  auto separated_volume_moments =
      cell_grouped_data_m->calculateNormalizedVolumeMoments(a_reconstruction);
  for (UnsignedIndex_t n = 0; n < 3; ++n) {
    guess_values_m(n) = separated_volume_moments[0].centroid()[n];
    guess_values_m(3 + n) = separated_volume_moments[1].centroid()[n];
  }
  guess_values_m = guess_values_m.cwiseProduct(weights_m);
}

template <class CellType>
const PlanarSeparator& MOFCommon<CellType>::getBestReconstruction(void) {
  return best_reconstruction_m;
}

template <class CellType>
const PlanarSeparator& MOFCommon<CellType>::getGuessReconstruction(void) {
  return guess_reconstruction_m;
}

template <class CellType>
double MOFCommon<CellType>::getInternalVolume(void) {
  return cell_grouped_data_m->getStoredMoments()[0].volume();
}

template <class CellType>
const Pt& MOFCommon<CellType>::getInternalCentroid(void) {
  return cell_grouped_data_m->getStoredMoments()[0].centroid();
}

template <class CellType>
const Pt& MOFCommon<CellType>::getExternalCentroid(void) {
  return cell_grouped_data_m->getStoredMoments()[1].centroid();
}

template <class CellType>
void MOF_2D<CellType>::setup(
    const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_data) {
  // Default setting of equal weighting for internal/external centroid.
  this->setup(a_cell_grouped_data, 0.5, 0.5);
}

template <class CellType>
void MOF_2D<CellType>::setup(
    const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_data,
    const double a_liquid_weight, const double a_gas_weight) {
  this->cell_grouped_data_m = &a_cell_grouped_data;
  // Need to do this way to prevent round off error from giving non-zero
  // z-component
  Normal centroid_line_normal = Normal::normalized(
      this->getExternalCentroid()[0] - this->getInternalCentroid()[0],
      this->getExternalCentroid()[1] - this->getInternalCentroid()[1], 0.0);
  if (magnitude(centroid_line_normal) < 0.9) {
    // Centroids provided a bad guess, start with a set normal
    centroid_line_normal =
        Normal(0.5 * std::sqrt(2.0), 0.5 * std::sqrt(2.0), 0.0);
  }
  Normal tangent0 = Normal(0.0, 0.0, 1.0);
  Normal tangent1 = crossProductNormalized(tangent0, centroid_line_normal);
  this->best_reference_frame_m =
      ReferenceFrame(tangent0, tangent1, centroid_line_normal);
  this->fillGeometryAndWeightVectors(a_liquid_weight, a_gas_weight);
}

template <class CellType>
PlanarSeparator MOF_2D<CellType>::solve(
    const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_moments,
    const double a_internal_weight, const double a_external_weight) {
  return this->runOptimization(this, a_cell_grouped_moments, a_internal_weight,
                               a_external_weight);
}

template <class CellType>
auto MOF_2D<CellType>::getDefaultInitialDelta(void)
    -> const Eigen::Matrix<double, columns_m, 1>& {
  initial_delta_m(0) = MOFCommon<CellType>::initial_angle_m;
  return initial_delta_m;
}

template <class CellType>
void MOF_2D<CellType>::updateGuess(
    Eigen::Matrix<double, columns_m, 1>* a_delta) {
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * (this->best_reference_frame_m);
  this->guess_reconstruction_m =
      this->getReconstructionFromParam(this->guess_reference_frame_m);
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
}

template <class CellType>
bool MOF_2D<CellType>::minimumReached(
    const Eigen::Matrix<double, columns_m, 1> a_delta) {
  return std::fabs(a_delta(0)) < MOFCommon<CellType>::minimum_angle_change_m;
}

template <class CellType>
inline UnitQuaternion MOF_2D<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

template <class CellType>
void MOF_3D<CellType>::setup(
    const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_data) {
  // Default setting of equal weighting for internal/external centroid.
  this->setup(a_cell_grouped_data, 0.5, 0.5);
}

template <class CellType>
void MOF_3D<CellType>::setup(
    const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_data,
    const double a_liquid_weight, const double a_gas_weight) {
  this->cell_grouped_data_m = &a_cell_grouped_data;
  Normal centroid_line_normal = Normal::fromPtNormalized(
      this->getExternalCentroid() - this->getInternalCentroid());
  if (magnitude(centroid_line_normal) < 0.9) {
    // Centroids provided a bad guess, start with a set normal
    centroid_line_normal = Normal(1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0),
                                  1.0 / std::sqrt(3.0));
  }
  this->best_reference_frame_m = getOrthonormalSystem(centroid_line_normal);
  this->fillGeometryAndWeightVectors(a_liquid_weight, a_gas_weight);
}

template <class CellType>
PlanarSeparator MOF_3D<CellType>::solve(
    const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_moments,
    const double a_internal_weight, const double a_external_weight) {
  return this->runOptimization(this, a_cell_grouped_moments, a_internal_weight,
                               a_external_weight);
}

template <class CellType>
auto MOF_3D<CellType>::getDefaultInitialDelta(void)
    -> const Eigen::Matrix<double, columns_m, 1>& {
  initial_delta_m(0) = MOFCommon<CellType>::initial_angle_m;
  initial_delta_m(1) = MOFCommon<CellType>::initial_angle_m;
  return initial_delta_m;
}

template <class CellType>
void MOF_3D<CellType>::updateGuess(
    Eigen::Matrix<double, columns_m, 1>* a_delta) {
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * (this->best_reference_frame_m);
  this->guess_reconstruction_m =
      this->getReconstructionFromParam(this->guess_reference_frame_m);
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
}

template <class CellType>
bool MOF_3D<CellType>::minimumReached(
    const Eigen::Matrix<double, columns_m, 1> a_delta) {
  return (std::fabs(a_delta(0)) < MOFCommon<CellType>::minimum_angle_change_m &&
          std::fabs(a_delta(1)) < MOFCommon<CellType>::minimum_angle_change_m);
}

template <class CellType>
inline UnitQuaternion MOF_3D<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(1), a_reference_frame[1]) *
         UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

template <class MOFType>
PlanarSeparator MOFDebug<MOFType>::solve(
    const CellGroupedMoments<typename MOFType::cell_type,
                             SeparatedMoments<VolumeMoments>>&
        a_cell_grouped_moments,
    const double a_internal_weight, const double a_external_weight) {
  return this->runOptimization(this, a_cell_grouped_moments, a_internal_weight,
                               a_external_weight);
}

template <class MOFType>
void MOFDebug<MOFType>::updateGuess(
    Eigen::Matrix<double, MOFType::columns_m, 1>* a_delta) {
  MOFType::updateGuess(a_delta);
  guess_reconstruction_history.push_back(MOFType::getGuessReconstruction());
}

template <class MOFType>
void MOFDebug<MOFType>::updateBestGuess(void) {
  MOFType::updateBestGuess();
  best_reconstruction_history.push_back(MOFType::getBestReconstruction());
}

template <class MOFType>
PlanarSeparator MOFDebug<MOFType>::getFinalReconstruction(void) {
  this->writeOutCentroidsAndWeights();
  for (std::size_t n = 0; n < best_reconstruction_history.size(); ++n) {
    this->writeOutPlane(best_reconstruction_history[n], "bestPlane", n);
  }
  return MOFType::getFinalReconstruction();
}

template <class MOFType>
void MOFDebug<MOFType>::writeOutPlane(const PlanarSeparator& a_reconstruction,
                                      const std::string& a_prefix,
                                      const std::size_t a_iteration_number) {
  for (UnsignedIndex_t plane = 0; plane < a_reconstruction.getNumberOfPlanes();
       ++plane) {
    Polygon poly = getPlanePolygonFromReconstruction<Polygon>(
        this->cell_grouped_data_m->getCell(), a_reconstruction,
        a_reconstruction[plane]);
    for (UnsignedIndex_t n = 0; n < poly.getNumberOfVertices(); ++n) {
      std::cout << a_prefix << "(1:3," << n + 1 << "," << 1 << "," << plane + 1
                << "," << a_iteration_number + 1 << ") = [";
      std::cout << poly[n].x() << "," << poly[n].y() << "," << poly[n].z()
                << "];\n";
    }
    std::cout << a_prefix << "Nvert(" << 1 << "," << plane + 1 << ","
              << a_iteration_number + 1 << ") = ";
    std::cout << poly.getNumberOfVertices() << ";\n";
  }
  std::cout << a_prefix << "NPlane(" << a_iteration_number + 1 << "," << 1
            << ") = ";
  std::cout << a_reconstruction.getNumberOfPlanes() << ";\n";
}

template <class MOFType>
void MOFDebug<MOFType>::writeOutCentroidsAndWeights(void) {
  for (int n = 0; n < 1; ++n) {
    std::cout << "liq_centroid(1:4," << n + 1 << ") = [";
    Pt liq_centroid(MOFType::correct_values_m(0), MOFType::correct_values_m(1),
                    MOFType::correct_values_m(2));
    liq_centroid = liq_centroid / safelyEpsilon(MOFType::weights_m(0));
    std::cout << liq_centroid.x() << "," << liq_centroid.y() << ","
              << liq_centroid.z() << "," << MOFType::weights_m(0) << "];\n";
    std::cout << "gas_centroid(1:4," << n + 1 << ") = [";
    Pt gas_centroid(MOFType::correct_values_m(3), MOFType::correct_values_m(4),
                    MOFType::correct_values_m(5));
    gas_centroid = gas_centroid / safelyEpsilon(MOFType::weights_m(3));
    std::cout << gas_centroid.x() << "," << gas_centroid.y() << ","
              << gas_centroid.z() << "," << MOFType::weights_m(3) << "];\n";
  }
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_MOF_TPP_
