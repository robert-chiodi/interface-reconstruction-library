// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_MOF_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_MOF_H_

#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>  // Eigen header

#include "src/generic_cutting/cut_polygon.h"
#include "src/geometry/general/rotations.h"
#include "src/geometry/polygons/polygon.h"
#include "src/interface_reconstruction_methods/volume_fraction_matching.h"
#include "src/moments/cell_grouped_moments.h"
#include "src/moments/separated_volume_moments.h"
#include "src/optimization/bfgs.h"
#include "src/optimization/levenberg_marquardt.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

// Forward declare to friend
template <class CellType>
class MOF_2D;

template <class CellType>
class MOF_3D;

template <class CellType>
class MOFCommon {
  friend MOF_2D<CellType>;
  friend MOF_3D<CellType>;

 public:
  /// \brief If `this->calculateScalarError()` is less than this, exit.
  static constexpr double acceptable_error_m = 1.0e-4 * 1.0e-4;
  /// \brief Maximum number of attempted iterations before exiting.
  static constexpr UnsignedIndex_t maximum_iterations_m = 40;
  /// \brief Minimum change in angle related delta below which minimum is
  /// deemed reached.
  static constexpr double minimum_angle_change_m = 0.0001 * 0.0174533;
  /// \brief Increase factor for lambda if more damping needed.
  static constexpr double lambda_increase_m = 5.0;
  /// \brief Decrease factor for lambda if new best solution is found.
  static constexpr double lambda_decrease_m = 1.0 / 10.0;
  /// \brief Number of iterations to allow between calculating a new Jacobian.
  static constexpr UnsignedIndex_t delay_jacobian_amount_m = 0;
  /// \brief Initial angle to use when first calculating Jacobian, equal to
  /// 5 degrees in radians.
  static constexpr double initial_angle_m =
      0.0001 * 0.0174533;  // 1e-4 Deg in radians

  /// \brief Default constructor
  MOFCommon(void) = default;

  template <class MOFType>
  PlanarSeparator runOptimization(
      MOFType* a_ptr_to_MOF_object,
      const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_moments,
      const double a_internal_weight, const double a_external_weight);

  /// \brief Return the final reconstruction to be used, which includes
  /// re-rotation and flipping of the reconstruction (while
  /// `getBestReconstruction()` and `getGuessReconstruction()` does not.
  PlanarSeparator getFinalReconstruction(void);

  /// \brief Calculate the error ||a_correct_vector - a_attempt_vector||^2
  /// where both vectors are weighted geometry vectors.
  double calculateScalarError(void);

  /// \brief Calculate the vector error a_correct_vector - a_attempt_vector
  /// where both vectors are weighted geometry vectors.
  Eigen::Matrix<double, 6, 1> calculateVectorError(void);

  /// \brief Calculate and return vector of differences between `guess_values_m`
  /// and `best_values_m`. Used in Jacobian calculation.
  Eigen::Matrix<double, 6, 1> calculateChangeInGuess(void);
  /// \brief Increase lambda to decrease step size
  void increaseLambda(double* a_lambda);

  /// \brief Decrease lambda to increase step size
  void decreaseLambda(double* a_lambda);

  /// \brief Return a boolean stating whether the value is too high or not.
  bool errorTooHigh(const double a_error);

  /// \brief Return a boolean stating whether maximum iterations has been
  /// exceeded.
  bool iterationTooHigh(const UnsignedIndex_t a_iteration);

  /// \brief Return whether Jacobian should be computed, allowing delayed
  /// re-evalulation of the Jacobian.
  bool shouldComputeJacobian(const UnsignedIndex_t a_iteration,
                             const UnsignedIndex_t a_last_jacobian);

  /// \brief Save current guess as the best guess
  void updateBestGuess(void);

  /// \brief Default destructor
  ~MOFCommon(void) = default;

 private:
  /// \brief Setup the weight and correct vectors.
  inline void fillGeometryAndWeightVectors(double a_liquid_weight,
                                           double a_gas_weight);

  /// \brief Perform geometric integration for volume moments, apply
  /// weighting, and save to guess_values_m.
  inline void setWeightedGeometryVectorFromReconstruction(
      const PlanarSeparator& a_reconstruction);

  /// \brief Given a reference frame, return the PlanarSeparator
  /// that matches the volume fraction with the normal from the reference frame.
  inline PlanarSeparator getReconstructionFromParam(
      const ReferenceFrame& a_reference_frame);

  /// \brief Return the best reconstruction found during the optimization
  /// procedure.
  const PlanarSeparator& getBestReconstruction(void);

  /// \brief Return a const reference to the guess reconstruction
  const PlanarSeparator& getGuessReconstruction(void);

  /// \brief Return the internal (under plane) volume from the
  /// cell_grouped_data_m.
  double getInternalVolume(void);

  /// \brief Return the internal (under plane) centroid from the
  /// cell_grouped_data_m.
  const Pt& getInternalCentroid(void);

  /// \brief Return the external (above plane) centroid from the
  /// cell_grouped_data_m.
  const Pt& getExternalCentroid(void);

  //------------------ Constants set during setup -------------------------
  const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>*
      cell_grouped_data_m;
  /// \brief Scale factor used for geometry in order to give dimension
  /// independence. Equal to vol^(1/3)
  double geom_scale_factor_m;
  /// \brief Weights to be applied to the corect and guess
  /// SeparatedMoments<VolumeMoments>/Surface Area.
  Eigen::Matrix<double, 6, 1> weights_m;
  /// \brief Weigted vector of correct values we are trying to match.
  Eigen::Matrix<double, 6, 1> correct_values_m;
  /// \brief Volume fraction in cell being reconstructed.
  double volume_fraction_m;
  //----------------------------------------------------------------------

  //---------------------- Working variables  ----------------------------
  /// \brief Weighted vector of guess values that we are trying to
  ///  match to `correct_values_m`.
  Eigen::Matrix<double, 6, 1> guess_values_m;
  /// \brief Guess reconstruction that is used to obtain `guess_values_m`.
  PlanarSeparator guess_reconstruction_m;
  /// \brief Guess reference frame used when obtaining `guess_reconstruction_m`.
  ReferenceFrame guess_reference_frame_m;
  //----------------------------------------------------------------------

  // Values saved throughout solution
  //-------------- Values saved throughout solution  ---------------------
  /// \brief Weighted vector that resulted in lowest error.
  Eigen::Matrix<double, 6, 1> best_values_m;
  /// \brief Reconstruction that has resulted in lowest error.
  PlanarSeparator best_reconstruction_m;
  /// \brief Reference frame associated with the best reconstruction.
  ReferenceFrame best_reference_frame_m;
  //----------------------------------------------------------------------
};

template <class CellType>
class MOF_2D : public MOFCommon<CellType> {
 public:
  using cell_type = CellType;
  /// \brief Number of rows for MOF_2D optimization.
  static constexpr int rows_m = 6;
  /// \brief Number of columns for MOF_2D optimization.
  static constexpr int columns_m = 1;

  /// \brief Default constructor
  MOF_2D(void) = default;

  void setup(
      const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_data);

  void setup(
      const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_data,
      const double a_liquid_weight, const double a_gas_weight);

  /// \brief Setup and solve the system, returning the found PlanarSeparator.
  PlanarSeparator solve(
      const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_moments,
      const double a_internal_weight, const double a_external_weight);

  auto getDefaultInitialDelta(void)
      -> const Eigen::Matrix<double, columns_m, 1>&;

  void updateGuess(Eigen::Matrix<double, columns_m, 1>* a_delta);

  bool minimumReached(const Eigen::Matrix<double, columns_m, 1> a_delta);

  /// \brief Default destructor
  ~MOF_2D(void) = default;

 private:
  inline UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);

  Eigen::Matrix<double, columns_m, 1> initial_delta_m;
};

template <class CellType>
class MOF_3D : public MOFCommon<CellType> {
 public:
  using cell_type = CellType;
  /// \brief Number of rows for MOF_2D optimization.
  static constexpr int rows_m = 6;
  /// \brief Number of columns for MOF_2D optimization.
  static constexpr int columns_m = 2;

  /// \brief Default constructor
  MOF_3D(void) = default;

  void setup(
      const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_data);

  void setup(
      const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_data,
      const double a_liquid_weight, const double a_gas_weight);

  /// \brief Setup and solve the system, returning the found PlanarSeparator.
  PlanarSeparator solve(
      const CellGroupedMoments<CellType, SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_moments,
      const double a_internal_weight, const double a_external_weight);

  auto getDefaultInitialDelta(void)
      -> const Eigen::Matrix<double, columns_m, 1>&;

  void updateGuess(Eigen::Matrix<double, columns_m, 1>* a_delta);

  bool minimumReached(const Eigen::Matrix<double, columns_m, 1> a_delta);

  /// \brief Default destructor
  ~MOF_3D(void) = default;

 private:
  inline UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);

  Eigen::Matrix<double, columns_m, 1> initial_delta_m;
};

/// \brief This class just calls the MOFType functions
/// but allows debug statements to be printed.
///
/// This class masks MOFType functions to allow calling
/// the function and then all printing any debug information
/// wanted. It also stores the guessed reconstructions and
/// saved "best" reconstructions in order to reprint them
/// out later. The printing format of this is to plot
/// the resulting vertices of the plane polygons
/// to be plotted by `R2P/references/r2p_history_plot.m`.
/// To do this, the information printed to screen needs to be
/// copied to `R2P/references/reconstruction_history.m`.
///
/// Requirements for MOFType class:
/// - MOFType class must meet all requirements
/// needed by the `LevenbergMarquardt` class in `optimizers.h`.
/// Right now, this is mainly MOF_2D and MOF_3d.
template <class MOFType>
class MOFDebug : public MOFType {
 public:
  /// \brief Default constructor
  MOFDebug(void) = default;

  /// \brief Setup and solve the system, returning the found PlanarSeparator.
  PlanarSeparator solve(
      const CellGroupedMoments<typename MOFType::cell_type,
                               SeparatedMoments<VolumeMoments>>&
          a_cell_grouped_moments,
      const double a_internal_weight, const double a_external_weight);

  /// \brief Calls updateGuess from the base class
  /// and stores the current reconstruction.
  void updateGuess(Eigen::Matrix<double, MOFType::columns_m, 1>* a_delta);

  /// \brief Calls updateBestGuess from the base class
  /// and stores the current reconstruction.
  void updateBestGuess(void);

  /// \brief Shadowed getFinalReconstruction call that writes out
  /// stored best reconstructions and then returns the reconstruction.
  PlanarSeparator getFinalReconstruction(void);

  /// \brief Default destructor
  ~MOFDebug(void) = default;

 private:
  /// \brief Write the ConvexPolygons in the reconstruction out
  /// to std::cout, tagged with the given iteration number.
  void writeOutPlane(const PlanarSeparator& a_reconstruction,
                     const std::string& a_prefix,
                     const std::size_t a_iteration_number);

  /// \brief Write out the centroids and weights to
  /// to enable visualization of what optimization is driving towards.
  void writeOutCentroidsAndWeights(void);

  /// \brief Saved guess reconstructions encountered during optimization.
  std::vector<PlanarSeparator> guess_reconstruction_history;
  /// \brief Saved best reconstructions accepted during optimization.
  std::vector<PlanarSeparator> best_reconstruction_history;
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/mof.tpp"
#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_MOF_H_
