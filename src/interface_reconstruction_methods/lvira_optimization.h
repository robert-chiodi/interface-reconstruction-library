// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_OPTIMIZATION_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_OPTIMIZATION_H_

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>  // Eigen header

#include "src/generic_cutting/cut_polygon.h"
#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/general/rotations.h"
#include "src/geometry/general/unit_quaternion.h"
#include "src/geometry/polygons/polygon.h"
#include "src/interface_reconstruction_methods/lvira_neighborhood.h"
#include "src/interface_reconstruction_methods/reconstruction_cleaning.h"
#include "src/optimization/bfgs.h"
#include "src/optimization/levenberg_marquardt.h"
#include "src/parameters/compiler_type.h"
#include "src/parameters/defined_types.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

/// \file lvira_optimization.h
///
/// This file contains all LVIRA functions needed
/// to be used during calls to the templated
/// Levenberg Marquardt solver in `optimizers.h`.

/// \brief Number of columns for 2D LVIRA optimization,
/// which is equal to parameters being fit.
static constexpr UnsignedIndex_t LVIRA_2D_columns = 1;

/// \brief Number of columns for 3D LVIRA optimization,
/// which is equal to parameters being fit.
static constexpr UnsignedIndex_t LVIRA_3D_columns = 2;

// Forward declaration to allow friending in LVIRACommon
template <class CellType>
class LVIRA_2D;
template <class CellType>
class LVIRA_3D;

/// \brief Class to contain data and methods that will be used
/// in all of the specific LVIRA optimization classes.
///
/// This class holds common methods and members that are used in the
/// specific LVIRA optimization classes (2D or 3D).
/// This includes variables and members that determine exit criterion
/// for the optimization through Levenberg-Marquardt.
///
/// Template parameters:
/// - `CellType` : Type of cell invovled in the neighborhood
/// - `kColumns` : Number of columns in Jacobian matrix. Equal to number of
/// parameters to optimize for with Levenberg-Marquardt.
template <class CellType, UnsignedIndex_t kColumns>
class LVIRACommon {
  template <class LVIRATypeForDebug>
  friend class LVIRADebug;
  friend LVIRA_2D<CellType>;
  friend LVIRA_3D<CellType>;

 public:
  /// \brief If `this->calculateScalarError()` is less than this, exit.
  static constexpr double acceptable_error_m = 1.0e-4 * 1.0e-4;
  /// \brief Maximum number of attempted iterations before exiting.
  static constexpr UnsignedIndex_t maximum_iterations_m = 20;
  /// \brief Minimum change in angle related delta below which minimum is
  /// deemed reached.
  static constexpr double minimum_angle_change_m = 0.0001745329;
  /// \brief Increase factor for lambda if more damping needed.
  static constexpr double lambda_increase_m = 5.0;
  /// \brief Decrease factor for lambda if new best solution is found.
  static constexpr double lambda_decrease_m = 1.0 / 10.0;
  /// \brief Number of iterations to allow between calculating a new Jacobian.
  static constexpr UnsignedIndex_t delay_jacobian_amount_m = 0;
  /// \brief Angle change to use when calculating finite-difference Jacobian.
  static constexpr double fininite_difference_angle_m =
      0.001 * 0.0174533;  // 1e-3 Deg in radians

  // Default constructor
  LVIRACommon(void) = default;

  /// \brief Class that takes a pointer to a LVIRA object
  /// and executes the optimization, returning the found
  /// optimal reconstruction.
  ///
  /// Requirements for LVIRAType:
  /// - Must be either `LVIRA_2D` or `LVIRA_3D`.
  ///
  /// \param[in] a_ptr_to_LVIRA_object Pointer to LVIRA object that will
  /// be run through optimization.
  /// \param[in] a_neighborhood_geometry Neighborhood geometry that will
  /// be used in the optimization to calculate error we are trying to minimize.
  /// \param[in] a_reconstruction Initial reconstruction to start LVIRA
  /// optimization routine from.
  template <class LVIRAType>
  PlanarSeparator runOptimization(
      LVIRAType* a_ptr_to_LVIRA_object,
      const LVIRANeighborhood<CellType>& a_neighborhood_geometry,
      const PlanarSeparator& a_reconstruction);

  /// \brief Return the final reconstruction to be used.
  PlanarSeparator getFinalReconstruction(void);

  /// \brief Calculate the vector error correct_values_m - guess_values_m
  /// where both vectors already have weight applied.
  Eigen::Matrix<double, Eigen::Dynamic, 1> calculateVectorError(void);

  /// \brief Return a boolean stating whether the value is too high or not.
  bool errorTooHigh(const double a_error);

  /// \brief Return a boolean stating whether maximum iterations has been
  /// exceeded.
  bool iterationTooHigh(const UnsignedIndex_t a_iteration);

  /// \brief Increase lambda to decrease step size
  void increaseLambda(double* a_lambda);

  /// \brief Decrease lambda to increase step size
  void decreaseLambda(double* a_lambda);

  /// \brief Return whether Jacobian should be computed, allowing delayed
  /// re-evalulation of the Jacobian.
  bool shouldComputeJacobian(const UnsignedIndex_t a_iteration,
                             const UnsignedIndex_t a_last_jacobian);

#ifndef USING_INTEL_COMPILER
  /// \brief Returns whether the minimum is reached.
  bool minimumReached(const Eigen::Matrix<double, static_cast<int>(kColumns),
                                          1>& a_delta) const;

  /// \brief Return vector of steps in parameters to take when
  /// computing Jacobian.
  Eigen::Matrix<double, static_cast<int>(kColumns), 1> getJacobianStepSize(
      void) const;
#else  // Intel really doesn't like the static cast of kColumns
  /// \brief Returns whether the minimum is reached.
  bool minimumReached(const Eigen::Matrix<double, kColumns, 1>& a_delta) const;

  /// \brief Return vector of steps in parameters to take when
  /// computing Jacobian.
  Eigen::Matrix<double, kColumns, 1> getJacobianStepSize(void) const;
#endif

  /// \brief Calculate the error ||correct_values_m - guess_values_m||^2
  /// where both vectors have already had weighting applied.
  double calculateScalarError(void);

  /// \brief Save current guess as the best guess.
  void updateBestGuess(void);

  /// \brief Calculate and return vector of differences between `guess_values_m`
  /// and `best_values_m`. Used in Jacobian calculation.
  Eigen::Matrix<double, Eigen::Dynamic, 1> calculateChangeInGuess(void);

  /// \brief Return the best reconstruction found during the optimization
  /// procedure.
  const PlanarSeparator& getBestReconstruction(void);

  /// \brief Return a const reference to the guess reconstruction
  const PlanarSeparator& getGuessReconstruction(void);

  /// \brief Return best guess reference frame
  const ReferenceFrame& getBestReferenceFrame(void);

  /// \brief Return the guess reference frame
  const ReferenceFrame& getGuessReferenceFrame(void);

  /// \brief Allocate the dynamically sized matrices.
  void allocateMatrices(const UnsignedIndex_t a_neighborhood_size);

  /// \brief Set up system and vectors with correct values/weights.
  void fillGeometryAndWeightVectors(void);

  /// \brief Update guess_values_m given a reconstruction.
  void setWeightedGeometryVectorFromReconstruction(
      const PlanarSeparator& a_reconstruction);

  /// \brief Return the corresponding planar reconstruction for a given
  /// reference frame.
  PlanarSeparator getReconstructionFromLVIRAParam(
      const ReferenceFrame& a_reference_frame);

  // Default destructor
  ~LVIRACommon(void) = default;
  // private: // TODO figure out how to make this private.

  // TODO Fix this and turn back to private. For some compilers, the LVIRA
  // functions complain they can't access the private members. Seen happen with
  // GNU 7.x.x compiler versions. Doesn't recognize the friend declaration?
 public:
  //------------------ Constants set during setup -------------------------
  /// \brief Pointer to neighborhood being used in reconstruction.
  const LVIRANeighborhood<CellType>* neighborhood_m;
  /// \brief Weights to be applied to the corect and guess values
  Eigen::Matrix<double, Eigen::Dynamic, 1> weights_m;
  /// \brief Weigted vector of correct values we are trying to match.
  Eigen::Matrix<double, Eigen::Dynamic, 1> correct_values_m;
  //----------------------------------------------------------------------

  //---------------------- Working variables  ----------------------------
  /// \brief Weighted vector of guess values that we are trying to
  ///  match to `correct_values_m`.
  Eigen::Matrix<double, Eigen::Dynamic, 1> guess_values_m;
  /// \brief Guess reconstruction that is used to obtain `guess_values_m`.
  PlanarSeparator guess_reconstruction_m;
  /// \brief Guess reference frame used when obtaining `guess_reconstruction_m`.
  ReferenceFrame guess_reference_frame_m;
  //----------------------------------------------------------------------

  // Values saved throughout solution
  //-------------- Values saved throughout solution  ---------------------
  /// \brief Weighted vector that resulted in lowest error.
  Eigen::Matrix<double, Eigen::Dynamic, 1> best_values_m;
  /// \brief Reconstruction that has resulted in lowest error.
  PlanarSeparator best_reconstruction_m;
  /// \brief Reference frame associated with the best reconstruction.
  ReferenceFrame best_reference_frame_m;
  //----------------------------------------------------------------------
};

/// \brief LVIRA class for reconstructions in 2D (x-y plane).
template <class CellType>
class LVIRA_2D : public LVIRACommon<CellType, LVIRA_2D_columns> {
 public:
  using cell_type = CellType;
  /// \brief Number of columns for LVIRA_2D optimization.
  static constexpr UnsignedIndex_t columns_m = LVIRA_2D_columns;

  /// \brief Default constructor
  LVIRA_2D(void) = default;
  /// \brief Default destructor
  ~LVIRA_2D(void) = default;

  /// \brief Perform optimization and return planar separator.
  PlanarSeparator solve(const LVIRANeighborhood<CellType>& a_neighborhood,
                        const PlanarSeparator& a_reconstruction);

  /// \brief Initialize simulation parameters that will
  /// be needed during the optimization.
  void setup(const PlanarSeparator& a_reconstruction);

  /// \brief Updates current guess reconstruction using a_delta
  /// and calculates the new reconstruction as well as
  /// stores weighted guess_value vector in `guess_values_m`.
  void updateGuess(const Eigen::Matrix<double, columns_m, 1>* const a_delta);

 private:
  /// \brief Return rotation for LVIRA dictated by elements in `a_delta`.
  ///
  /// The rotation order is:
  /// - Rotate by `a_delta(0)` radians around `a_reference_frame[0]`
  UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);
};

/// \brief LVIRA class for reconstructions in 3 dimensions.
template <class CellType>
class LVIRA_3D : public LVIRACommon<CellType, LVIRA_3D_columns> {
 public:
  using cell_type = CellType;
  /// \brief Number of columns for LVIRA_3D optimization.
  static constexpr UnsignedIndex_t columns_m = LVIRA_3D_columns;

  /// \brief Default constructor
  LVIRA_3D(void) = default;
  /// \brief Default destructor
  ~LVIRA_3D(void) = default;

  /// \brief Perform optimization and return planar separator.
  PlanarSeparator solve(const LVIRANeighborhood<CellType>& a_neighborhood,
                        const PlanarSeparator& a_reconstruction);

  /// \brief Initialize simulation parameters that will
  /// be needed during the optimization.
  void setup(const PlanarSeparator& a_reconstruction);

  /// \brief Updates current guess reconstruction using a_delta
  /// and calculates the new reconstruction as well as
  /// stores weighted guess_value vector in `guess_values_m`.
  void updateGuess(const Eigen::Matrix<double, columns_m, 1>* const a_delta);

 private:
  /// \brief Return rotation for LVIRA dictated by elements in `a_delta`.
  ///
  /// The rotation order is:
  /// - Rotate by `a_delta(0)` radians around `a_reference_frame[0]`
  /// - Rotate by `a_delta(1)` radians around `a_reference_frame[1]`
  UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);
};

/// \brief This class just calls the LVIRAType functions
/// but allows debug statements to be printed. The
/// solution path is also saved to be exported and
/// visualized.
template <class LVIRAType>
class LVIRADebug : public LVIRAType {
 public:
  /// \brief Default constructor
  LVIRADebug(void) = default;

  /// \brief Calls updateGuess from the base class
  /// and stores the current reconstruction.
  void updateGuess(
      const Eigen::Matrix<double, LVIRAType::columns_m, 1>* const a_delta);

  /// \brief Solve call for debugging LVIRA.
  PlanarSeparator solve(
      const LVIRANeighborhood<typename LVIRAType::cell_type>& a_neighborhood,
      const PlanarSeparator& a_reconstruction);

  /// \brief Calls updateBestGuess from the base class
  /// and stores the current best reconstruction.
  void updateBestGuess(void);

  /// \brief Shadowed getFinalReconstruction call that writes out
  /// stored best reconstructions and then returns the reconstruction.
  PlanarSeparator getFinalReconstruction(void);

  /// \brief Write the Polygons in the reconstruction out
  /// to std::cout, tagged with the given iteration number.
  void writeOutPlane(const PlanarSeparator& a_reconstruction,
                     const std::string& a_prefix,
                     const std::size_t a_iteration_number);

  /// \brief Write out the volume fractions and weights to
  /// to enable visualization of what optimization is driving towards.
  void writeOutVolumesAndWeights(void);

  /// \brief Default destructor
  ~LVIRADebug(void) = default;

 private:
  /// \brief Saved guess reconstructions encountered during optimization.
  /// Note: This will also have reconstructions caused by the
  /// calculation of the Jacobian.
  std::vector<PlanarSeparator> guess_reconstruction_history;
  /// \brief Saved best reconstructions accepted during optimization.
  std::vector<PlanarSeparator> best_reconstruction_history;
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/lvira_optimization.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_OPTIMIZATION_H_
