// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_OPTIMIZATION_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_OPTIMIZATION_TPP_

namespace IRL {

template <class CellType, UnsignedIndex_t kColumns>
template <class R2PType>
PlanarSeparator R2PCommon<CellType, kColumns>::runOptimization(
    R2PType* a_ptr_to_R2P_object,
    const R2PNeighborhood<CellType>& a_neighborhood_geometry,
    const PlanarSeparator& a_reconstruction) {
  a_ptr_to_R2P_object->setup(a_neighborhood_geometry, a_reconstruction);
  LevenbergMarquardt<R2PType, -1, static_cast<int>(R2PType::columns_m)>
      lm_solver;
  lm_solver.solve(a_ptr_to_R2P_object,
                  static_cast<int>(correct_values_m.rows()),
                  a_ptr_to_R2P_object->getDefaultInitialDelta());
  //  BFGS<R2PType, static_cast<int>(R2PType::columns_m)> bfgs_solver;
  //  bfgs_solver.solve(a_ptr_to_R2P_object,
  //                    a_ptr_to_R2P_object->getDefaultInitialDelta());
  return a_ptr_to_R2P_object->getFinalReconstruction();
}

template <class CellType, UnsignedIndex_t kColumns>
PlanarSeparator R2PCommon<CellType, kColumns>::getFinalReconstruction(void) {
  PlanarSeparator reconstruction = this->getBestReconstruction();
  cleanReconstruction(system_center_cell_m, volume_fraction_m, &reconstruction);
  for (auto& plane : reconstruction) {
    plane.distance() += plane.normal() * initial_center_cell_centroid_m;
  }
  return reconstruction;
}

template <class CellType, UnsignedIndex_t kColumns>
Eigen::Matrix<double, Eigen::Dynamic, 1>
R2PCommon<CellType, kColumns>::calculateVectorError(void) {
  return (correct_values_m - guess_values_m);
}

template <class CellType, UnsignedIndex_t kColumns>
bool R2PCommon<CellType, kColumns>::errorTooHigh(const double a_error) {
  return a_error > acceptable_error_m;
}

template <class CellType, UnsignedIndex_t kColumns>
bool R2PCommon<CellType, kColumns>::iterationTooHigh(
    const UnsignedIndex_t a_iteration) {
  return a_iteration > maximum_iterations_m;
}

template <class CellType, UnsignedIndex_t kColumns>
void R2PCommon<CellType, kColumns>::increaseLambda(double* a_lambda) {
  (*a_lambda) *= lambda_increase_m;
}

template <class CellType, UnsignedIndex_t kColumns>
void R2PCommon<CellType, kColumns>::decreaseLambda(double* a_lambda) {
  (*a_lambda) *= lambda_decrease_m;
}

template <class CellType, UnsignedIndex_t kColumns>
bool R2PCommon<CellType, kColumns>::shouldComputeJacobian(
    const UnsignedIndex_t a_iteration, const UnsignedIndex_t a_last_jacobian) {
  return a_iteration - a_last_jacobian > delay_jacobian_amount_m;
}

template <class CellType, UnsignedIndex_t kColumns>
double R2PCommon<CellType, kColumns>::calculateScalarError(void) {
  return (correct_values_m - guess_values_m).squaredNorm();
}

template <class CellType, UnsignedIndex_t kColumns>
void R2PCommon<CellType, kColumns>::saveBestGuess(void) {
  best_values_m = guess_values_m;
  best_reconstruction_m = guess_reconstruction_m;
  best_reference_frame_m = guess_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kColumns>
Eigen::Matrix<double, Eigen::Dynamic, 1>
R2PCommon<CellType, kColumns>::calculateChangeInGuess(void) {
  return guess_values_m - best_values_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const PlanarSeparator& R2PCommon<CellType, kColumns>::getBestReconstruction(
    void) {
  return best_reconstruction_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const PlanarSeparator& R2PCommon<CellType, kColumns>::getGuessReconstruction(
    void) {
  return guess_reconstruction_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const CellType& R2PCommon<CellType, kColumns>::getCell(
    const UnsignedIndex_t a_cell) {
  return cells_to_cut_m[a_cell];
}

template <class CellType, UnsignedIndex_t kColumns>
const ReferenceFrame& R2PCommon<CellType, kColumns>::getBestReferenceFrame(
    void) {
  return best_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kColumns>
const ReferenceFrame& R2PCommon<CellType, kColumns>::getGuessReferenceFrame(
    void) {
  return guess_reference_frame_m;
}

template <class CellType, UnsignedIndex_t kColumns>
void R2PCommon<CellType, kColumns>::addCellMomentsToGeometryVector(
    const SeparatedMoments<VolumeMoments>& a_moments_from_cut_cell,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& a_weights,
    const UnsignedIndex_t a_start_index,
    Eigen::Matrix<double, Eigen::Dynamic, 1>* a_geometry_vector) {
  (*a_geometry_vector)(a_start_index) =
      a_weights(a_start_index) * a_moments_from_cut_cell[0].volume();
  (*a_geometry_vector)(a_start_index + 1) =
      a_weights(a_start_index + 1) * a_moments_from_cut_cell[0].centroid()[0];
  (*a_geometry_vector)(a_start_index + 2) =
      a_weights(a_start_index + 2) * a_moments_from_cut_cell[0].centroid()[1];
  (*a_geometry_vector)(a_start_index + 3) =
      a_weights(a_start_index + 3) * a_moments_from_cut_cell[0].centroid()[2];
  (*a_geometry_vector)(a_start_index + 4) =
      a_weights(a_start_index + 4) * a_moments_from_cut_cell[1].centroid()[0];
  (*a_geometry_vector)(a_start_index + 5) =
      a_weights(a_start_index + 5) * a_moments_from_cut_cell[1].centroid()[1];
  (*a_geometry_vector)(a_start_index + 6) =
      a_weights(a_start_index + 6) * a_moments_from_cut_cell[1].centroid()[2];
}

template <class CellType, UnsignedIndex_t kColumns>
void R2PCommon<CellType, kColumns>::setWeightedGeometryVectorFromReconstruction(
    const PlanarSeparator& a_reconstruction) {
  for (UnsignedIndex_t i = 0; i < cells_to_cut_m.size(); ++i) {
    UnsignedIndex_t start_index = 7 * i;
    auto moments_from_cut_cell =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>,
                                   ReconstructionDefaultCuttingMethod>(
            cells_to_cut_m[i], a_reconstruction);
    addCellMomentsToGeometryVector(moments_from_cut_cell, weights_m,
                                   start_index, &guess_values_m);
  }
}

template <class CellType, UnsignedIndex_t kColumns>
inline void
R2PCommon<CellType, kColumns>::setWeightedGeometryVectorFromSurfaceArea(
    const PlanarSeparator& a_reconstruction) {
  guess_values_m(guess_values_m.rows() - 1) =
      weights_m(guess_values_m.rows() - 1) *
      std::sqrt(
          getReconstructionSurfaceArea(system_center_cell_m, a_reconstruction));
}

// Turn off warnings about sign conversion because need to work
// with Eigen which using long int, and vector which uses std::size_t
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
template <class CellType, UnsignedIndex_t kColumns>
void R2PCommon<CellType, kColumns>::fillGeometryAndWeightVectors(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const double distance_multiplier, const double volume_weight_switch) {
  // Get and store correct center centroids
  const Pt& center_liquid_centroid =
      (a_neighborhood.getCenterCellStoredMoments())[0].centroid();
  const Pt& center_gas_centroid =
      (a_neighborhood.getCenterCellStoredMoments())[1].centroid();

  // Fill cells_to_cut_m, correct_values_m, and weights_m.
  double liquid_volume_weight_magnitude = 0.0;
  double liquid_centroid_weight_magnitude = 0.0;
  double gas_centroid_weight_magnitude = 0.0;
  stencil_average_volume_fraction_m = 0.0;

  // Save the initial center cell centroid before shifting everyting to center
  // around 0.0.
  initial_center_cell_centroid_m =
      a_neighborhood.getCenterCell().calculateCentroid();

  double average_cell_volume = 0.0;
  cells_to_cut_m.resize(a_neighborhood.size());
  for (UnsignedIndex_t n = 0; n < a_neighborhood.size(); ++n) {
    // Add cell
    cells_to_cut_m[n] = a_neighborhood.getCell(n);
    for (UnsignedIndex_t v = 0; v < cells_to_cut_m[n].getNumberOfVertices();
         ++v) {
      cells_to_cut_m[n][v] -= initial_center_cell_centroid_m;
    }
    average_cell_volume += cells_to_cut_m[n].calculateVolume();
  }
  system_center_cell_m = a_neighborhood.getCenterCell();
  for (UnsignedIndex_t v = 0; v < system_center_cell_m.getNumberOfVertices();
       ++v) {
    system_center_cell_m[v] -= initial_center_cell_centroid_m;
  }

  average_cell_volume /= static_cast<double>(a_neighborhood.size());
  characteristic_length_m = std::pow(average_cell_volume, 1.0 / 3.0);

  for (UnsignedIndex_t n = 0; n < a_neighborhood.size(); ++n) {
    const Pt& liquid_centroid =
        (a_neighborhood.getStoredMoments(n))[0].centroid();
    const Pt& gas_centroid = (a_neighborhood.getStoredMoments(n))[1].centroid();

    double liquid_centroid_distance2 =
        squaredDistanceBetweenPts(liquid_centroid, center_liquid_centroid) /
        (characteristic_length_m * characteristic_length_m);
    double gas_centroid_distance2 =
        squaredDistanceBetweenPts(gas_centroid, center_gas_centroid) /
        (characteristic_length_m * characteristic_length_m);

    // Add correct moments
    correct_values_m(7 * n + 0) =
        a_neighborhood.getStoredMoments(n)[0].volume();
    correct_values_m(7 * n + 1) =
        liquid_centroid[0] - initial_center_cell_centroid_m[0];
    correct_values_m(7 * n + 2) =
        liquid_centroid[1] - initial_center_cell_centroid_m[1];
    correct_values_m(7 * n + 3) =
        liquid_centroid[2] - initial_center_cell_centroid_m[2];
    correct_values_m(7 * n + 4) =
        gas_centroid[0] - initial_center_cell_centroid_m[0];
    correct_values_m(7 * n + 5) =
        gas_centroid[1] - initial_center_cell_centroid_m[1];
    correct_values_m(7 * n + 6) =
        gas_centroid[2] - initial_center_cell_centroid_m[2];

    // Add even weighting of 1 for liquid volume
    weights_m(7 * n) = 1.0;

    double cell_volume = cells_to_cut_m[n].calculateVolume();

    // Add weighting for liquid centroid if there is liquid
    if (correct_values_m(7 * n) / cell_volume > global_constants::VF_LOW) {
      double front_multiplier = {(1.0 - volume_weight_switch) +
                                 correct_values_m(7 * n) / cell_volume *
                                     volume_weight_switch};
      weights_m(7 * n + 1) =
          front_multiplier *
          std::exp(-distance_multiplier * liquid_centroid_distance2);
    } else {
      weights_m(7 * n + 1) = 0.0;
    }
    weights_m(7 * n + 2) = weights_m(7 * n + 1);
    weights_m(7 * n + 3) = weights_m(7 * n + 1);
    // Add weighting for gas centroid if there is gas
    if (correct_values_m(7 * n) / cell_volume < global_constants::VF_HIGH) {
      double front_multiplier = {(1.0 - volume_weight_switch) +
                                 (1.0 - correct_values_m(7 * n) / cell_volume) *
                                     volume_weight_switch};
      weights_m(7 * n + 4) =
          front_multiplier *
          std::exp(-distance_multiplier * gas_centroid_distance2);
    } else {
      weights_m(7 * n + 4) = 0.0;
    }
    weights_m(7 * n + 5) = weights_m(7 * n + 4);
    weights_m(7 * n + 6) = weights_m(7 * n + 4);

    // Track average volume fraction in stencil
    stencil_average_volume_fraction_m +=
        (correct_values_m(7 * n) / cell_volume);

    // Keep track of the magnitudes of each kind of weight
    liquid_volume_weight_magnitude += weights_m(7 * n) * weights_m(7 * n);
    liquid_centroid_weight_magnitude +=
        weights_m(7 * n + 1) * weights_m(7 * n + 1);
    gas_centroid_weight_magnitude +=
        weights_m(7 * n + 4) * weights_m(7 * n + 4);
  }

  // Get average volume fraction
  stencil_average_volume_fraction_m /=
      static_cast<double>(a_neighborhood.size());

  // Normalize each kind of weight to have a L2 magnitude of 1 now.
  liquid_volume_weight_magnitude =
      1.0 / std::sqrt(liquid_volume_weight_magnitude);
  liquid_centroid_weight_magnitude =
      1.0 / std::sqrt(liquid_centroid_weight_magnitude);
  gas_centroid_weight_magnitude =
      1.0 / std::sqrt(gas_centroid_weight_magnitude);

  for (UnsignedIndex_t n = 0; n < guess_values_m.rows() - 1; n += 7) {
    weights_m(n) *= liquid_volume_weight_magnitude;
    weights_m(n + 1) *= liquid_centroid_weight_magnitude;
    weights_m(n + 2) *= liquid_centroid_weight_magnitude;
    weights_m(n + 3) *= liquid_centroid_weight_magnitude;
    weights_m(n + 4) *= gas_centroid_weight_magnitude;
    weights_m(n + 5) *= gas_centroid_weight_magnitude;
    weights_m(n + 6) *= gas_centroid_weight_magnitude;
  }

  // Add surface area with weighting of 1.0.
  correct_values_m(guess_values_m.rows() - 1) =
      std::sqrt(a_neighborhood.getSurfaceArea());
  weights_m(guess_values_m.rows() - 1) = 1.0;
}
#pragma GCC diagnostic pop

template <class CellType, UnsignedIndex_t kColumns>
void R2PCommon<CellType, kColumns>::setRelativeImportanceBetweenWeights(
    double a_importance_of_liquid_volume_fraction,
    double a_importance_of_liquid_centroid_relative_to_gas,
    double a_importance_of_centroid, double a_importance_of_surface_area) {
  assert(a_importance_of_liquid_volume_fraction >= 0.0);
  assert(a_importance_of_liquid_centroid_relative_to_gas >= 0.0 &&
         a_importance_of_liquid_centroid_relative_to_gas <= 1.0);
  assert(a_importance_of_centroid >= 0.0);
  assert(a_importance_of_surface_area >= 0.0);

  double sum = a_importance_of_liquid_volume_fraction +
               a_importance_of_centroid + a_importance_of_surface_area;

  double importance_of_liquid_centroid = {
      a_importance_of_liquid_centroid_relative_to_gas *
      a_importance_of_centroid};
  double importance_of_gas_centroid = {
      (1.0 - a_importance_of_liquid_centroid_relative_to_gas) *
      a_importance_of_centroid};

  // Make weights sum to 1.
  a_importance_of_liquid_volume_fraction /= sum;
  importance_of_liquid_centroid /= sum;
  importance_of_gas_centroid /= sum;
  a_importance_of_surface_area /= sum;
  assert(a_importance_of_liquid_volume_fraction +
             importance_of_liquid_centroid + importance_of_gas_centroid +
             a_importance_of_surface_area >
         1.0 - 1.0e-13);
  // Add on scaling factor
  a_importance_of_liquid_volume_fraction /=
      (characteristic_length_m * characteristic_length_m *
       characteristic_length_m);
  importance_of_liquid_centroid /= characteristic_length_m;
  importance_of_gas_centroid /= characteristic_length_m;
  a_importance_of_surface_area /= characteristic_length_m;

  // Readjust weights
  for (UnsignedIndex_t n = 0; n < weights_m.rows() - 1; n += 7) {
    weights_m(n) *= a_importance_of_liquid_volume_fraction;
    weights_m(n + 1) *= importance_of_liquid_centroid;
    weights_m(n + 2) *= importance_of_liquid_centroid;
    weights_m(n + 3) *= importance_of_liquid_centroid;
    weights_m(n + 4) *= importance_of_gas_centroid;
    weights_m(n + 5) *= importance_of_gas_centroid;
    weights_m(n + 6) *= importance_of_gas_centroid;
  }
  weights_m(weights_m.rows() - 1) *= a_importance_of_surface_area;
}

template <class R2PType>
void R2PDebug<R2PType>::updateGuess(
    Eigen::Matrix<double, R2PType::columns_m, 1>* a_delta) {
  R2PType::updateGuess(a_delta);
  guess_reconstruction_history.push_back(R2PType::getGuessReconstruction());
}

template <class R2PType>
PlanarSeparator R2PDebug<R2PType>::solve(
    const R2PNeighborhood<typename R2PType::cell_type>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class R2PType>
void R2PDebug<R2PType>::updateBestGuess(void) {
  R2PType::updateBestGuess();
  best_reconstruction_history.push_back(R2PType::getBestReconstruction());
  PlanarSeparator tmp = best_reconstruction_history.back();
  std::cout << "**************************************** \n";
  std::cout << "Writing out plane: " << best_reconstruction_history.size() - 1
            << '\n';
  std::cout << tmp << std::endl;
}

template <class R2PType>
PlanarSeparator R2PDebug<R2PType>::getFinalReconstruction(void) {
  this->writeOutCentroidsAndWeights();
  for (std::size_t n = 0; n < best_reconstruction_history.size(); ++n) {
    this->writeOutPlane(best_reconstruction_history[n], "bestPlane", n);
  }
  return R2PType::getFinalReconstruction();
}

template <class R2PType>
void R2PDebug<R2PType>::writeOutPlane(const PlanarSeparator& a_reconstruction,
                                      const std::string& a_prefix,
                                      const std::size_t a_iteration_number) {
  for (UnsignedIndex_t cell = 0; cell < this->cells_to_cut_m.size(); ++cell) {
    for (UnsignedIndex_t plane = 0;
         plane < a_reconstruction.getNumberOfPlanes(); ++plane) {
      Polygon poly = getPlanePolygonFromReconstruction<Polygon>(
          this->cells_to_cut_m[cell], a_reconstruction,
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

template <class R2PType>
void R2PDebug<R2PType>::writeOutCentroidsAndWeights(void) {
  for (UnsignedIndex_t n = 0; n < this->cells_to_cut_m.size(); ++n) {
    std::cout << "liq_centroid(1:4," << n + 1 << ") = [";
    Pt liq_centroid(R2PType::correct_values_m(7 * n + 1),
                    R2PType::correct_values_m(7 * n + 2),
                    R2PType::correct_values_m(7 * n + 3));
    liq_centroid = liq_centroid / safelyEpsilon(R2PType::weights_m(7 * n + 1));
    // liq_centroid += (-1.0 * this->initial_center_cell_centroid_m);
    std::cout << liq_centroid.x() << "," << liq_centroid.y() << ","
              << liq_centroid.z() << "," << R2PType::weights_m(7 * n + 1)
              << "];\n";
    std::cout << "gas_centroid(1:4," << n + 1 << ") = [";
    Pt gas_centroid(R2PType::correct_values_m(7 * n + 4),
                    R2PType::correct_values_m(7 * n + 5),
                    R2PType::correct_values_m(7 * n + 6));
    gas_centroid = gas_centroid / safelyEpsilon(R2PType::weights_m(7 * n + 4));
    // gas_centroid += (-1.0 * this->initial_center_cell_centroid_m);
    std::cout << gas_centroid.x() << "," << gas_centroid.y() << ","
              << gas_centroid.z() << "," << R2PType::weights_m(7 * n + 4)
              << "];\n";
  }
}

//******************************************************************* //
//******************************************************************* //
//******************************************************************* //
// Individual R2P functions below this
//******************************************************************* //
//******************************************************************* //
//******************************************************************* //

//******************************************************************* //
//     Functions for R2P_2D1P below this
//******************************************************************* //
/// \brief Perform optimization and return planar separator.
template <class CellType>
PlanarSeparator R2P_2D1P<CellType>::solve(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class CellType>
void R2P_2D1P<CellType>::setup(const R2PNeighborhood<CellType>& a_neighborhood,
                               const PlanarSeparator& a_reconstruction) {
  this->allocateMatrices(a_neighborhood.size());
  this->volume_fraction_m =
      a_neighborhood.getCenterCellStoredMoments()[0].volume() /
      a_neighborhood.getCenterCell().calculateVolume();

  double distance_multiplier = {1.0};
  double volume_weight_switch = {1.0};
  this->fillGeometryAndWeightVectors(a_neighborhood, distance_multiplier,
                                     volume_weight_switch);
  this->setRelativeImportanceBetweenWeights(0.0, 0.5, 1.0, 0.0);
  this->correct_values_m = this->correct_values_m.cwiseProduct(this->weights_m);

  Normal tangent0 = Normal(0.0, 0.0, 1.0);
  Normal tangent1 =
      crossProductNormalized(tangent0, a_reconstruction[0].normal());
  this->best_reference_frame_m =
      ReferenceFrame(tangent0, tangent1, a_reconstruction[0].normal());
}

template <class CellType>
void R2P_2D1P<CellType>::allocateMatrices(
    const UnsignedIndex_t a_neighborhood_size) {
  const int rows = static_cast<int>(a_neighborhood_size * 7) + 1;
  this->weights_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->correct_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->guess_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->best_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
}

template <class CellType>
void R2P_2D1P<CellType>::updateGuess(
    const Eigen::Matrix<double, columns_m, 1>* const a_delta) {
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * this->best_reference_frame_m;
  this->guess_reconstruction_m =
      this->getReconstructionFromR2PParam(this->guess_reference_frame_m);
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
  this->setWeightedGeometryVectorFromSurfaceArea(this->guess_reconstruction_m);
}

template <class CellType>
void R2P_2D1P<CellType>::updateBestGuess(void) {
  this->saveBestGuess();
}

template <class CellType>
UnitQuaternion R2P_2D1P<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

template <class CellType>
PlanarSeparator R2P_2D1P<CellType>::getReconstructionFromR2PParam(
    const ReferenceFrame& a_reference_frame) {
  // Distance guess is 0.0 because that is a plane going through center of cell.
  // The cell was specifically moved to have cell centroid at (0.0, 0.0, 0.0).
  auto planar_separator =
      PlanarSeparator::fromOnePlane(Plane(a_reference_frame[2], 0.0));
  setDistanceToMatchVolumeFractionPartialFill(
      this->system_center_cell_m, this->volume_fraction_m, &planar_separator);
  return planar_separator;
}

/// \brief Returns whether the minimum is reached.
template <class CellType>
bool R2P_2D1P<CellType>::minimumReached(
    const Eigen::Matrix<double, columns_m, 1>& delta) {
  return std::fabs(delta(0)) <
         R2PCommon<CellType, columns_m>::minimum_angle_change_m;
}

/// \brief Return initial delta to be used in optimization
template <class CellType>
auto R2P_2D1P<CellType>::getDefaultInitialDelta(void)
    -> const Eigen::Matrix<double, columns_m, 1>& {
  initial_delta_m(0) = R2PCommon<CellType, columns_m>::initial_angle_m;
  return initial_delta_m;
}

//******************************************************************* //
//     Functions for R2P_3D1P below this
//******************************************************************* //
template <class CellType>
PlanarSeparator R2P_3D1P<CellType>::solve(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class CellType>
void R2P_3D1P<CellType>::setup(const R2PNeighborhood<CellType>& a_neighborhood,
                               const PlanarSeparator& a_reconstruction) {
  this->allocateMatrices(a_neighborhood.size());
  this->volume_fraction_m =
      a_neighborhood.getCenterCellStoredMoments()[0].volume() /
      a_neighborhood.getCenterCell().calculateVolume();

  double distance_multiplier = {1.0};
  double volume_weight_switch = {1.0};
  this->fillGeometryAndWeightVectors(a_neighborhood, distance_multiplier,
                                     volume_weight_switch);
  this->setRelativeImportanceBetweenWeights(0.0, 0.5, 1.0, 0.0);
  this->correct_values_m = this->correct_values_m.cwiseProduct(this->weights_m);

  this->best_reference_frame_m =
      getOrthonormalSystem(a_reconstruction[0].normal());
}

template <class CellType>
void R2P_3D1P<CellType>::allocateMatrices(
    const UnsignedIndex_t a_neighborhood_size) {
  const int rows = static_cast<int>(a_neighborhood_size * 7) + 1;
  this->weights_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->correct_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->guess_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->best_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
}

template <class CellType>
void R2P_3D1P<CellType>::updateGuess(
    const Eigen::Matrix<double, columns_m, 1>* const a_delta) {
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * this->best_reference_frame_m;
  this->guess_reconstruction_m =
      this->getReconstructionFromR2PParam(this->guess_reference_frame_m);
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
  this->setWeightedGeometryVectorFromSurfaceArea(this->guess_reconstruction_m);
}

template <class CellType>
void R2P_3D1P<CellType>::updateBestGuess(void) {
  this->saveBestGuess();
}

template <class CellType>
bool R2P_3D1P<CellType>::minimumReached(
    const Eigen::Matrix<double, columns_m, 1>& delta) {
  return std::max(std::fabs(delta(0)), std::fabs(delta(1))) <
         R2PCommon<CellType, columns_m>::minimum_angle_change_m;
}

/// \brief Return initial delta to be used in optimization
template <class CellType>
auto R2P_3D1P<CellType>::getDefaultInitialDelta(void)
    -> const Eigen::Matrix<double, columns_m, 1>& {
  initial_delta_m(0) = R2PCommon<CellType, columns_m>::initial_angle_m;
  initial_delta_m(1) = R2PCommon<CellType, columns_m>::initial_angle_m;
  return initial_delta_m;
}

template <class CellType>
UnitQuaternion R2P_3D1P<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(1), a_reference_frame[1]) *
         UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

template <class CellType>
PlanarSeparator R2P_3D1P<CellType>::getReconstructionFromR2PParam(
    const ReferenceFrame& a_reference_frame) {
  // Distance guess is 0.0 because that is a plane going through center of cell.
  // The cell was specifically moved to have cell centroid at (0.0, 0.0, 0.0).
  auto planar_separator =
      PlanarSeparator::fromOnePlane(Plane(a_reference_frame[2], 0.0));
  setDistanceToMatchVolumeFractionPartialFill(
      this->system_center_cell_m, this->volume_fraction_m, &planar_separator);
  return planar_separator;
}

//******************************************************************* //
//     Functions for R2P_2D2P below this
//******************************************************************* //
template <class CellType>
PlanarSeparator R2P_2D2P<CellType>::solve(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class CellType>
void R2P_2D2P<CellType>::setup(const R2PNeighborhood<CellType>& a_neighborhood,
                               const PlanarSeparator& a_reconstruction) {
  this->allocateMatrices(a_neighborhood.size());
  this->volume_fraction_m =
      a_neighborhood.getCenterCellStoredMoments()[0].volume() /
      a_neighborhood.getCenterCell().calculateVolume();

  double distance_multiplier = {1.0};
  double volume_weight_switch = {1.0};
  this->fillGeometryAndWeightVectors(a_neighborhood, distance_multiplier,
                                     volume_weight_switch);
  double weight_surface_area =
      0.75 *
      std::pow(2.0 * (0.5 - this->stencil_average_volume_fraction_m), 2.0);
  this->setRelativeImportanceBetweenWeights(0.0, 0.5, 1.0, weight_surface_area);
  this->correct_values_m = this->correct_values_m.cwiseProduct(this->weights_m);

  Normal tangent0;
  Normal shared_normal =
      getSharedNormal(a_reconstruction[1].normal(),
                      a_reconstruction[0].normal(), &best_beta_m, &tangent0);

  tangent0 = Normal(0.0, 0.0, 1.0);  // Hard set since known in 2D

  // Calculate correct beta according to being flipped
  if (a_reconstruction.isFlipped()) {
    best_beta_m = angleNormalize(-best_beta_m);
  }

  Normal tangent1 = crossProductNormalized(tangent0, shared_normal);
  this->best_reference_frame_m =
      ReferenceFrame(tangent0, tangent1, shared_normal);
  Normal plus_beta_normal =
      UnitQuaternion(best_beta_m, tangent0) * shared_normal;
  this->best_distances_m[0] =
      a_reconstruction[0].distance() -
      a_reconstruction[0].normal() * this->initial_center_cell_centroid_m;
  this->best_distances_m[1] =
      a_reconstruction[1].distance() -
      a_reconstruction[1].normal() * this->initial_center_cell_centroid_m;
  if (plus_beta_normal != a_reconstruction[0].normal()) {
    std::swap(best_distances_m[0], best_distances_m[1]);
  }
}

template <class CellType>
void R2P_2D2P<CellType>::allocateMatrices(
    const UnsignedIndex_t a_neighborhood_size) {
  const int rows = static_cast<int>(a_neighborhood_size * 7 + 1);
  this->weights_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->correct_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->guess_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->best_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
}

template <class CellType>
void R2P_2D2P<CellType>::updateGuess(
    Eigen::Matrix<double, columns_m, 1>* a_delta) {
  // Rotate around z-axis stored in reference_frame_m[1]
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * this->best_reference_frame_m;
  this->guess_distances_m[0] = this->best_distances_m[0] + (*a_delta)(2);
  this->guess_distances_m[1] = this->best_distances_m[1] + (*a_delta)(3);
  this->guess_beta_m = angleNormalize(best_beta_m + (*a_delta)(1));
  this->guess_reconstruction_m = this->getReconstructionFromR2PParam(
      this->guess_reference_frame_m, guess_beta_m, guess_distances_m);
  guess_distances_m[0] = this->guess_reconstruction_m[0].distance();
  guess_distances_m[1] = this->guess_reconstruction_m[1].distance();
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
  this->setWeightedGeometryVectorFromSurfaceArea(this->guess_reconstruction_m);
  // Adjust to have correct gradient when calculating Jacobian
  (*a_delta)(2) = (guess_distances_m[0] - best_distances_m[0]) /
                  this->characteristic_length_m;
  (*a_delta)(3) = (guess_distances_m[1] - best_distances_m[1]) /
                  this->characteristic_length_m;
}

template <class CellType>
void R2P_2D2P<CellType>::updateBestGuess(void) {
  this->saveBestGuess();
  best_beta_m = guess_beta_m;
  best_distances_m[0] = guess_distances_m[0];
  best_distances_m[1] = guess_distances_m[1];
}

template <class CellType>
bool R2P_2D2P<CellType>::minimumReached(
    const Eigen::Matrix<double, columns_m, 1>& delta) {
  return (std::fabs(delta(0)) <
              R2PCommon<CellType, columns_m>::minimum_angle_change_m &&
          std::fabs(delta(1)) <
              R2PCommon<CellType, columns_m>::minimum_angle_change_m &&
          std::fabs(delta(2)) <
              R2PCommon<CellType, columns_m>::minimum_distance_change_m *
                  this->characteristic_length_m &&
          std::fabs(delta(3)) <
              R2PCommon<CellType, columns_m>::minimum_distance_change_m *
                  this->characteristic_length_m);
}

template <class CellType>
auto R2P_2D2P<CellType>::getDefaultInitialDelta(void)
    -> const Eigen::Matrix<double, columns_m, 1>& {
  initial_delta_m(0) = R2PCommon<CellType, columns_m>::initial_angle_m;
  initial_delta_m(1) = R2PCommon<CellType, columns_m>::initial_angle_m;
  initial_delta_m(2) = R2PCommon<CellType, columns_m>::initial_distance_m *
                       this->characteristic_length_m;
  initial_delta_m(3) = R2PCommon<CellType, columns_m>::initial_distance_m *
                       this->characteristic_length_m;
  return initial_delta_m;
}

template <class CellType>
const double& R2P_2D2P<CellType>::getBestBeta(void) {
  return best_beta_m;
}

template <class CellType>
UnitQuaternion R2P_2D2P<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

template <class CellType>
PlanarSeparator R2P_2D2P<CellType>::getReconstructionFromR2PParam(
    const ReferenceFrame& a_reference_frame, const double a_beta,
    const double* a_distances) {
  UnitQuaternion beta_rotation(a_beta, a_reference_frame[0]);
  Plane plane_0(beta_rotation * a_reference_frame[2], a_distances[0]);
  Plane plane_1((beta_rotation.inverse()) * a_reference_frame[2],
                a_distances[1]);
  double iswitch = angleNormalize(a_beta) < M_PI ? 1.0 : -1.0;
  PlanarSeparator reconstruction_to_return =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, iswitch);
  setDistanceToMatchVolumeFractionPartialFill(this->system_center_cell_m,
                                              this->volume_fraction_m,
                                              &reconstruction_to_return);
  cleanReconstructionSameNormal(this->system_center_cell_m,
                                this->volume_fraction_m,
                                &reconstruction_to_return);
  return reconstruction_to_return;
}

//******************************************************************* //
//     Functions for R2P_3D2P below this
//******************************************************************* //
template <class CellType>
PlanarSeparator R2P_3D2P<CellType>::solve(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_reconstruction) {
  return this->runOptimization(this, a_neighborhood, a_reconstruction);
}

template <class CellType>
void R2P_3D2P<CellType>::setup(const R2PNeighborhood<CellType>& a_neighborhood,
                               const PlanarSeparator& a_reconstruction) {
  this->allocateMatrices(a_neighborhood.size());
  this->volume_fraction_m =
      a_neighborhood.getCenterCellStoredMoments()[0].volume() /
      a_neighborhood.getCenterCell().calculateVolume();

  double distance_multiplier = {1.0};
  double volume_weight_switch = {1.0};
  this->fillGeometryAndWeightVectors(a_neighborhood, distance_multiplier,
                                     volume_weight_switch);
  double weight_surface_area =
      0.75 *
      std::pow(2.0 * (0.5 - this->stencil_average_volume_fraction_m), 2.0);
  this->setRelativeImportanceBetweenWeights(0.0, 0.5, 1.0, weight_surface_area);
  this->correct_values_m = this->correct_values_m.cwiseProduct(this->weights_m);

  Normal tangent0;
  Normal shared_normal =
      getSharedNormal(a_reconstruction[1].normal(),
                      a_reconstruction[0].normal(), &best_beta_m, &tangent0);

  // If cutting is flipped, beta should be > pi by convention
  if (a_reconstruction.isFlipped()) {
    best_beta_m = angleNormalize(-best_beta_m);
  }

  Normal tangent1 = crossProductNormalized(tangent0, shared_normal);
  this->best_reference_frame_m =
      ReferenceFrame(tangent0, tangent1, shared_normal);
  Normal plus_beta_normal =
      UnitQuaternion(best_beta_m, tangent0) * shared_normal;
  best_distances_m[0] =
      a_reconstruction[0].distance() -
      a_reconstruction[0].normal() * this->initial_center_cell_centroid_m;
  best_distances_m[1] =
      a_reconstruction[1].distance() -
      a_reconstruction[1].normal() * this->initial_center_cell_centroid_m;
  if (plus_beta_normal != a_reconstruction[0].normal()) {
    std::swap(best_distances_m[0], best_distances_m[1]);
  }
}

template <class CellType>
void R2P_3D2P<CellType>::allocateMatrices(
    const UnsignedIndex_t a_neighborhood_size) {
  const int rows = static_cast<int>(a_neighborhood_size * 7 + 1);
  this->weights_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->correct_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->guess_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
  this->best_values_m = Eigen::Matrix<double, Eigen::Dynamic, 1>(rows, 1);
}

template <class CellType>
void R2P_3D2P<CellType>::updateGuess(
    Eigen::Matrix<double, columns_m, 1>* a_delta) {
  UnitQuaternion rotation =
      this->getDeltaRotationQuat(this->best_reference_frame_m, *a_delta);
  this->guess_reference_frame_m = rotation * this->best_reference_frame_m;
  guess_distances_m[0] = best_distances_m[0] + (*a_delta)(4);
  guess_distances_m[1] = best_distances_m[1] + (*a_delta)(5);

  guess_beta_m = angleNormalize(best_beta_m + (*a_delta)(3));
  this->guess_reconstruction_m = this->getReconstructionFromR2PParam(
      this->guess_reference_frame_m, guess_beta_m, guess_distances_m);
  guess_distances_m[0] = this->guess_reconstruction_m[0].distance();
  guess_distances_m[1] = this->guess_reconstruction_m[1].distance();
  this->setWeightedGeometryVectorFromReconstruction(
      this->guess_reconstruction_m);
  this->setWeightedGeometryVectorFromSurfaceArea(this->guess_reconstruction_m);
  // Adjust to have correct gradient when calculating Jacobian
  (*a_delta)(4) = (guess_distances_m[0] - best_distances_m[0]) /
                  this->characteristic_length_m;
  (*a_delta)(5) = (guess_distances_m[1] - best_distances_m[1]) /
                  this->characteristic_length_m;
}

template <class CellType>
void R2P_3D2P<CellType>::updateBestGuess(void) {
  this->saveBestGuess();
  best_beta_m = guess_beta_m;
  best_distances_m[0] = guess_distances_m[0];
  best_distances_m[1] = guess_distances_m[1];
}

template <class CellType>
bool R2P_3D2P<CellType>::minimumReached(
    const Eigen::Matrix<double, columns_m, 1>& delta) {
  return (std::fabs(delta(0)) <
              R2PCommon<CellType, columns_m>::minimum_angle_change_m &&
          std::fabs(delta(1)) <
              R2PCommon<CellType, columns_m>::minimum_angle_change_m &&
          std::fabs(delta(2)) <
              R2PCommon<CellType, columns_m>::minimum_angle_change_m &&
          std::fabs(delta(3)) <
              R2PCommon<CellType, columns_m>::minimum_angle_change_m &&
          std::fabs(delta(4)) <
              R2PCommon<CellType, columns_m>::minimum_distance_change_m *
                  this->characteristic_length_m &&
          std::fabs(delta(5)) <
              R2PCommon<CellType, columns_m>::minimum_distance_change_m *
                  this->characteristic_length_m);
}

template <class CellType>
auto R2P_3D2P<CellType>::getDefaultInitialDelta(void)
    -> const Eigen::Matrix<double, columns_m, 1>& {
  initial_delta_m(0) = R2PCommon<CellType, columns_m>::initial_angle_m;
  initial_delta_m(1) = R2PCommon<CellType, columns_m>::initial_angle_m;
  initial_delta_m(2) = R2PCommon<CellType, columns_m>::initial_angle_m;
  initial_delta_m(3) = R2PCommon<CellType, columns_m>::initial_angle_m;
  initial_delta_m(4) = R2PCommon<CellType, columns_m>::initial_distance_m *
                       this->characteristic_length_m;
  initial_delta_m(5) = R2PCommon<CellType, columns_m>::initial_distance_m *
                       this->characteristic_length_m;
  return initial_delta_m;
}

template <class CellType>
const double& R2P_3D2P<CellType>::getBestBeta(void) {
  return best_beta_m;
}

template <class CellType>
UnitQuaternion R2P_3D2P<CellType>::getDeltaRotationQuat(
    const ReferenceFrame& a_reference_frame,
    const Eigen::Matrix<double, columns_m, 1>& a_delta) {
  return UnitQuaternion(a_delta(2), a_reference_frame[2]) *
         UnitQuaternion(a_delta(1), a_reference_frame[1]) *
         UnitQuaternion(a_delta(0), a_reference_frame[0]);
}

template <class CellType>
PlanarSeparator R2P_3D2P<CellType>::getReconstructionFromR2PParam(
    const ReferenceFrame& a_reference_frame, const double a_beta,
    const double* a_distances) {
  UnitQuaternion beta_rotation(a_beta, a_reference_frame[0]);
  Plane plane_0(beta_rotation * a_reference_frame[2], a_distances[0]);
  Plane plane_1((beta_rotation.inverse()) * a_reference_frame[2],
                a_distances[1]);
  double iswitch = angleNormalize(a_beta) < M_PI ? 1.0 : -1.0;
  PlanarSeparator reconstruction_to_return =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, iswitch);
  setDistanceToMatchVolumeFractionPartialFill(this->system_center_cell_m,
                                              this->volume_fraction_m,
                                              &reconstruction_to_return);
  cleanReconstructionSameNormal(this->system_center_cell_m,
                                this->volume_fraction_m,
                                &reconstruction_to_return);
  return reconstruction_to_return;
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_OPTIMIZATION_TPP_
