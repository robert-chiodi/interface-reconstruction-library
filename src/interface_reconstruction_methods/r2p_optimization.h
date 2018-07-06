// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_OPTIMIZATION_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_OPTIMIZATION_H_

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
#include "src/geometry/polygons/polygon.h"
#include "src/interface_reconstruction_methods/r2p_neighborhood.h"
#include "src/interface_reconstruction_methods/reconstruction_cleaning.h"
#include "src/optimization/bfgs.h"
#include "src/optimization/levenberg_marquardt.h"
#include "src/parameters/defined_types.h"

namespace IRL {

/// \file r2p_optimization.h
///
/// This file contains all R2P functions needed
/// to be used during calls to the templated
/// Levenberg Marquardt solver in `optimizers.h`.
///
/// First, the function declarations are given.
/// Afterwards, the inlined function definitions are given.
/// Then, templated functions are given.

/// \brief Number of cells for 2D - 1 plane R2P optimization.
static constexpr UnsignedIndex_t R2P_2D1P_ncells = 9;
/// \brief Number of rows for 2D - 1 plane R2P optimization, with
/// 7 Entries per cell (Liquid volume, liquid centroid, gas centroid)
/// and the surface area for the cell being reconstructed.
static constexpr UnsignedIndex_t R2P_2D1P_rows = 7 * R2P_2D1P_ncells + 1;
/// \brief Number of columns for 2D - 1 plane R2P optimization,
/// which is equal to parameters being fit.
static constexpr UnsignedIndex_t R2P_2D1P_columns = 1;

/// \brief Number of cells for 3D - 1 plane R2P optimization.
static constexpr UnsignedIndex_t R2P_3D1P_ncells = 27;
/// \brief Number of rows for 3D - 1 plane R2P optimization, with
/// 7 Entries per cell (Liquid volume, liquid centroid, gas centroid)
/// and the surface area for the cell being reconstructed.
static constexpr UnsignedIndex_t R2P_3D1P_rows = 7 * R2P_3D1P_ncells + 1;
/// \brief Number of columns for 3D - 1 plane R2P optimization,
/// which is equal to parameters being fit.
static constexpr UnsignedIndex_t R2P_3D1P_columns = 2;

/// \brief Number of cells for 2D - 2 plane R2P optimization.
static constexpr UnsignedIndex_t R2P_2D2P_ncells = 9;
/// \brief Number of rows for 2D - 2 plane R2P optimization, with
/// 7 Entries per cell (Liquid volume, liquid centroid, gas centroid)
/// and the surface area for the cell being reconstructed.
static constexpr UnsignedIndex_t R2P_2D2P_rows = 7 * R2P_2D2P_ncells + 1;
/// \brief Number of columns for 2D - 2 plane R2P optimization,
/// which is equal to parameters being fit.
static constexpr UnsignedIndex_t R2P_2D2P_columns = 4;

/// \brief Number of cells for 3D - 2 plane R2P optimization.
static constexpr UnsignedIndex_t R2P_3D2P_ncells = 27;
/// \brief Number of rows for 3D - 2 plane R2P optimization, with
/// 7 Entries per cell (Liquid volume, liquid centroid, gas centroid)
/// and the surface area for the cell being reconstructed.
static constexpr UnsignedIndex_t R2P_3D2P_rows = 7 * R2P_3D2P_ncells + 1;
/// \brief Number of columns for 3D - 2 plane R2P optimization,
/// which is equal to parameters being fit.
static constexpr UnsignedIndex_t R2P_3D2P_columns = 6;

template <class CellType>
class R2P_2D1P;
template <class CellType>
class R2P_3D1P;
template <class CellType>
class R2P_2D2P;
template <class CellType>
class R2P_3D2P;

/// \brief Class to contain data and methods that will be used
/// in all of the specific R2Poptimization.
///
/// This class holds common methods and members that are used in the
/// specific R2P optimization classes (R2P_2D1P, R2P_3D1P, R2P_2D2P, R2P_3D2P).
/// This includes variables and members that determine exit criterion
/// for the optimization through Levenberg-Marquardt.
///
/// Template parameters:
/// - `kNcells` : Number of cells that will be included during the
/// optimization.
/// - `kRows` : Number of rows in the error/guess/correct/weight vectors.
/// - `kColumns` : Number of columns in Jacobian matrix. Equal to number of
/// parameters to optimize for with Levenberg-Marquardt.
template <class CellType, UnsignedIndex_t kColumns>
class R2PCommon {
  template <class R2PTypeForDebug>
  friend class R2PDebug;
  friend R2P_2D1P<CellType>;
  friend R2P_3D1P<CellType>;
  friend R2P_2D2P<CellType>;
  friend R2P_3D2P<CellType>;

 public:
  /// \brief If `this->calculateScalarError()` is less than this, exit.
  static constexpr double acceptable_error_m = 1.0e-4 * 1.0e-4;
  /// \brief Maximum number of attempted iterations before exiting.
  static constexpr UnsignedIndex_t maximum_iterations_m = 20;
  /// \brief Minimum change in angle related delta below which minimum is
  /// deemed reached.
  static constexpr double minimum_angle_change_m = 0.0001745329;
  /// \brief Minimum change in distance related delta below which minimum is
  /// deemed reached.
  static constexpr double minimum_distance_change_m = 1.0e-4;
  /// \brief Increase factor for lambda if more damping needed.
  static constexpr double lambda_increase_m = 5.0;
  /// \brief Decrease factor for lambda if new best solution is found.
  static constexpr double lambda_decrease_m = 1.0 / 10.0;
  /// \brief Number of iterations to allow between calculating a new Jacobian.
  static constexpr UnsignedIndex_t delay_jacobian_amount_m = 0;
  /// \brief Initial angle to use when first calculating Jacobian, equal to
  /// 5 degrees in radians.
  static constexpr double initial_angle_m =
      0.001 * 0.0174533;  // 1e-3 Deg in radians
  /// \brief Initial distance to use when first calculating Jacobian.
  static constexpr double initial_distance_m = 0.001;

  // Default constructor
  R2PCommon(void) = default;

  /// \brief Class that takes a pointer to a R2P object
  /// and executes the optimization, returning the found
  /// optimal reconstruction.
  ///
  /// Requirements for R2PType:
  /// - Must be either `R2P_2D1P`, `R2P_3D1P`, `R2P_2D2P`, or `R2P_3D2P`.
  ///
  /// \param[in] a_ptr_to_R2P_object Pointer to R2P object that will
  /// be run through optimization.
  /// \param[in] a_neighborhood_geometry Neighborhood geometry that will
  /// be used in the optimization to calculate error we are trying to minimize.
  /// \param[in] a_reconstruction Initial reconstruction to start R2P
  /// optimization routine from.
  template <class R2PType>
  PlanarSeparator runOptimization(
      R2PType* a_ptr_to_R2P_object,
      const R2PNeighborhood<CellType>& a_neighborhood_geometry,
      const PlanarSeparator& a_reconstruction);

  /// \brief Return the final reconstruction to be used.
  PlanarSeparator getFinalReconstruction(void);

  /// \brief Calculate the vector error a_correct_vector - a_attempt_vector
  /// where both vectors are weighted geometry vectors.
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

  /// \brief Calculate the error ||a_correct_vector - a_attempt_vector||^2 where
  /// both vectors are weighted geometry vectors.
  double calculateScalarError(void);

  /// \brief Save current `guess_values_m` and `current_reconstruction_m`.
  void saveBestGuess(void);

  /// \brief Calculate and return vector of differences between `guess_values_m`
  /// and `best_values_m`. Used in Jacobian calculation.
  Eigen::Matrix<double, Eigen::Dynamic, 1> calculateChangeInGuess(void);

  /// \brief Return the best reconstruction found during the optimization
  /// procedure.
  const PlanarSeparator& getBestReconstruction(void);

  /// \brief Return a const reference to the guess reconstruction
  const PlanarSeparator& getGuessReconstruction(void);

  /// \brief Return const reference to cells_to_cut
  const CellType& getCell(const UnsignedIndex_t a_cell);

  /// \brief Return best guess reference frame
  const ReferenceFrame& getBestReferenceFrame(void);

  /// \brief Return the guess reference frame
  const ReferenceFrame& getGuessReferenceFrame(void);

  // Default destructor
  ~R2PCommon(void) = default;
  // private: // TODO figure out how to make this private.

  /// \brief A convenience class that takes SeparatedMoments<VolumeMoments>
  /// after cutting a cell and fills geometry vector with it, starting at
  /// a_start_index.
  ///
  /// This function is to be used to place SeparatedMoments<VolumeMoments> from
  /// cut cells into a geometry vector to be used by R2P. Starting at
  /// `a_start_index`, it places 7 elements (Liquid volume, Liquid Centroid, Gas
  /// Centroid) into `a_geometry_vector`, which is a Eigen::Matrix vector. The
  /// template `ReturnGeometryVector` is needed to accept Eigen::Matrix vectors
  /// of different fixed sizes.
  ///
  /// Template Requirements for `ReturnGeometryVector`:
  /// - An accessing function `operator()` that accepts one argument
  /// (an integer as an index) and will save the incoming double upon
  /// assignment.
  ///
  /// \param[in] a_moments_from_cut_cell SeparatedMoments<VolumeMoments> to be
  /// added to a_geometry_vector \param[in] a_weights Vector of weights to
  /// component-wise multiply with `a_geometry_vector`. \param[inout]
  /// a_start_index Starting location to begin adding elements to
  /// `a_geometry_vector`. \param[out] a_geometry_vector Geometry vector to have
  /// elements added to.
  static void addCellMomentsToGeometryVector(
      const SeparatedMoments<VolumeMoments>& a_moments_from_cut_cell,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& a_weights,
      const UnsignedIndex_t a_start_index,
      Eigen::Matrix<double, Eigen::Dynamic, 1>* a_geometry_vector);

  /// \brief Initial filling of correct geometry and weighting vectors
  /// during setup.
  ///
  /// Assumes that Neighborhood geometry (all centroids/cells) have already
  /// been moved to be with respect to the center of the cell undergoing
  /// reconstruction (which would be at the 0,0,0 location in the stencil)
  /// and scaled so that each cell is a shifted unit cube. Upon returning,
  /// the correct_values_m vector will be filled with 7 entries per cell:
  /// - Liquid volume fraction
  /// - Liquid centroid x component
  /// - Liquid centroid y component
  /// - Liquid centroid z component
  /// - Gas centroid x component
  /// - Gas centroid y component
  /// - Gas centroid z component
  ///
  /// These entries will be in the "positive and descending" normal
  /// reference frame.
  void fillGeometryAndWeightVectors(
      const R2PNeighborhood<CellType>& a_neighborhood,
      const double distance_multiplier, const double volume_weight_switch);

  /// \brief Adjusts weight vectors so that the magnitude sums to one
  /// and satisfies the given magnitudes of the specific weights.
  ///
  /// \param[in] a_importance_of_liquid_volume_fraction
  /// How important liquid volume fraction is in the error.
  /// \param[in] a_importance_of_liquid_centroid_relative_to_gas
  /// How important liquid centroid is compared to the gas centroid
  /// \param[in] a_importance_of_centroid How important centroids are
  /// compared to liquid volume fraction and surface area.
  /// \param[in] a_importance_of_surface_area Importance of matching surface
  /// area.
  void setRelativeImportanceBetweenWeights(
      double a_importance_of_liquid_volume_fraction,
      double a_importance_of_liquid_centroid_relative_to_gas,
      double a_importance_of_centroid, double a_importance_of_surface_area);

  /// \brief Use reconstruction to calculate weighted `guess_values_m vector`
  void setWeightedGeometryVectorFromReconstruction(
      const PlanarSeparator& a_reconstruction);

  /// \brief Use reconstruction to calculate weighted surface area contribution
  /// to `guess_values_m`.
  inline void setWeightedGeometryVectorFromSurfaceArea(
      const PlanarSeparator& a_reconstruction);

  // TODO Fix this and turn back to private. For some compilers, the R2P
  // functions complain they can't access the private members. Seen happen with
  // GNU 7.x.x compiler versions.
 public:
  //------------------ Constants set during setup -------------------------
  /// \brief Cells involved in the optimization given by vertices.
  std::vector<CellType> cells_to_cut_m;
  /// \brief Weights to be applied to the corect and guess
  /// SeparatedMoments<VolumeMoments>/Surface Area.
  Eigen::Matrix<double, Eigen::Dynamic, 1> weights_m;
  /// \brief Weigted vector of correct values we are trying to match.
  Eigen::Matrix<double, Eigen::Dynamic, 1> correct_values_m;
  /// \brief Volume fraction in cell being reconstructed.
  double volume_fraction_m;
  /// \brief Average volume fraction in stencil, used to set some weightings.
  double stencil_average_volume_fraction_m;
  /// \brief Characteristic length to be used for normalization.
  double characteristic_length_m;
  /// \brief Center cell.
  CellType system_center_cell_m;
  /// \brief Center cell centroid to be used for shifting.
  Pt initial_center_cell_centroid_m;
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

/// \brief R2P class for reconstructions in 2 dimensions with 1 plane (hence
/// 2D1P)
template <class CellType>
class R2P_2D1P : public R2PCommon<CellType, R2P_2D1P_columns> {
 public:
  using cell_type = CellType;
  /// \brief Number of columns for R2P_2D1P optimization.
  static constexpr UnsignedIndex_t columns_m = R2P_2D1P_columns;

  /// \brief Default constructor
  R2P_2D1P(void) = default;
  /// \brief Default destructor
  ~R2P_2D1P(void) = default;

  /// \brief Perform optimization and return planar separator.
  PlanarSeparator solve(const R2PNeighborhood<CellType>& a_neighborhood,
                        const PlanarSeparator& a_reconstruction);

  /// \brief Initialize simulation parameters that will
  /// be needed during the optimization.
  void setup(const R2PNeighborhood<CellType>& a_neighborhood,
             const PlanarSeparator& a_reconstruction);

  void allocateMatrices(const UnsignedIndex_t a_neighborhood_size);

  /// \brief Updates current guess reconstruction using a_delta
  /// and calculates the new reconstruction as well as
  /// stores weighted guess_value vector in `guess_values_m`.
  ///
  /// Delta may need to be modified as well when placing the
  /// the planes requires reprojection of the distances
  /// onto a volume-conserving distance.
  void updateGuess(const Eigen::Matrix<double, columns_m, 1>* const a_delta);

  /// \brief Save current guess as best and update all internally
  /// saved variables such as reference frame, beta, etc.
  void updateBestGuess(void);

  /// \brief Returns whether the minimum is reached.
  bool minimumReached(const Eigen::Matrix<double, columns_m, 1>& delta);

  /// \brief Return initial delta to be used in optimization
  auto getDefaultInitialDelta(void)
      -> const Eigen::Matrix<double, columns_m, 1>&;

 private:
  /// \brief Return rotation for R2P dictated by elements in `a_delta`.
  ///
  /// The rotation order is:
  /// - Rotate by `a_delta(0)` radians around `a_reference_frame[0]`
  UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);

  /// \brief Return a reconstruction from R2P parameters
  PlanarSeparator getReconstructionFromR2PParam(
      const ReferenceFrame& a_reference_frame);

  /// \brief Initial delta to supply to Levenberg Marquardt
  ///  to calculate Jacobian.
  Eigen::Matrix<double, columns_m, 1> initial_delta_m;
};

/// \brief R2P class for reconstructions in 3 dimensions with 1 plane (hence
/// 3D1P)
template <class CellType>
class R2P_3D1P : public R2PCommon<CellType, R2P_3D1P_columns> {
 public:
  using cell_type = CellType;
  /// \brief Number of columns for R2P_3D1P optimization.
  static constexpr UnsignedIndex_t columns_m = R2P_3D1P_columns;

  /// \brief Default constructor
  R2P_3D1P(void) = default;
  /// \brief Default destructor
  ~R2P_3D1P(void) = default;

  /// \brief Perform optimization and return planar separator.
  PlanarSeparator solve(const R2PNeighborhood<CellType>& a_neighborhood,
                        const PlanarSeparator& a_reconstruction);

  /// \brief Initialize simulation parameters that will
  /// be needed during the optimization.
  void setup(const R2PNeighborhood<CellType>& a_neighborhood,
             const PlanarSeparator& a_reconstruction);

  void allocateMatrices(const UnsignedIndex_t a_neighborhood_size);

  /// \brief Updates current guess reconstruction using a_delta
  /// and calculates the new reconstruction as well as
  /// stores weighted guess_value vector in `guess_values_m`.
  ///
  /// Delta may need to be modified as well when placing the
  /// the planes requires reprojection of the distances
  /// onto a volume-conserving distance.
  void updateGuess(const Eigen::Matrix<double, columns_m, 1>* const a_delta);

  /// \brief Save current guess as best and update all internally
  /// saved variables such as reference frame, beta, etc.
  void updateBestGuess(void);

  /// \brief Returns whether the minimum is reached.
  bool minimumReached(const Eigen::Matrix<double, columns_m, 1>& delta);

  /// \brief Return initial delta to be used in optimization
  auto getDefaultInitialDelta(void)
      -> const Eigen::Matrix<double, columns_m, 1>&;

 private:
  /// \brief Return rotation for R2P dictated by elements in `a_delta`.
  ///
  /// The rotation order is:
  /// - Rotate by `a_delta(0)` radians around `a_reference_frame[0]`
  /// - Rotate by `a_delta(1)` radians around `a_reference_frame[1]`
  UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);

  /// \brief Return a reconstruction from R2P parameters
  PlanarSeparator getReconstructionFromR2PParam(
      const ReferenceFrame& a_reference_frame);

  /// \brief Initial delta to supply to Levenberg Marquardt
  ///  to calculate Jacobian.
  Eigen::Matrix<double, columns_m, 1> initial_delta_m;
};

/// \brief R2P class for reconstructions in 2 dimensions with 2 plane (hence
/// 2D2P)
template <class CellType>
class R2P_2D2P : public R2PCommon<CellType, R2P_2D2P_columns> {
 public:
  using cell_type = CellType;
  /// \brief Number of columns for R2P_2D2P optimization.
  static constexpr UnsignedIndex_t columns_m = R2P_2D2P_columns;

  /// \brief Default constructor
  R2P_2D2P(void) = default;
  /// \brief Default destructor
  ~R2P_2D2P(void) = default;

  /// \brief Perform optimization and return planar separator.
  PlanarSeparator solve(const R2PNeighborhood<CellType>& a_neighborhood,
                        const PlanarSeparator& a_reconstruction);

  /// \brief Initialize simulation parameters that will
  /// be needed during the optimization.
  void setup(const R2PNeighborhood<CellType>& a_neighborhood,
             const PlanarSeparator& a_reconstruction);

  void allocateMatrices(const UnsignedIndex_t a_neighborhood_size);

  /// \brief Updates current guess reconstruction using a_delta
  /// and calculates the new reconstruction as well as
  /// stores weighted guess_value vector in `guess_values_m`.
  ///
  /// Delta may need to be modified as well when placing the
  /// the planes requires reprojection of the distances
  /// onto a volume-conserving distance.
  void updateGuess(Eigen::Matrix<double, columns_m, 1>* a_delta);

  /// \brief Save current guess as best and update all internally
  /// saved variables such as reference frame, beta, etc.
  void updateBestGuess(void);

  /// \brief Returns whether the minimum is reached.
  bool minimumReached(const Eigen::Matrix<double, columns_m, 1>& delta);

  /// \brief Return initial delta to be used in optimization
  auto getDefaultInitialDelta(void)
      -> const Eigen::Matrix<double, columns_m, 1>&;

  /// \brief Return best beta_m value
  const double& getBestBeta(void);

 private:
  /// \brief Return rotation for R2P dictated by elements in `a_delta`.
  ///
  /// The rotation order is:
  /// - Rotate by `a_delta(0)` radians around `a_reference_frame[0]`
  UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);

  /// \brief Return a reconstruction from R2P parameters
  PlanarSeparator getReconstructionFromR2PParam(
      const ReferenceFrame& a_reference_frame, const double a_beta,
      const double* a_distances);

  /// \brief Initial delta to supply to Levenberg Marquardt
  ///  to calculate Jacobian.
  Eigen::Matrix<double, columns_m, 1> initial_delta_m;
  /// \brief Current guess for beta, the rotation from shared normal
  ///  to two plane normals.
  double guess_beta_m;
  /// \brief Current guess for distance to each plane.
  double guess_distances_m[2];
  /// \brief Current best value for beta, the rotation from shared normal
  ///  to two plane normals, that gives the least error.
  double best_beta_m;
  /// \brief Current best values for distance to each plane that gives
  /// the least error.
  double best_distances_m[2];
};

/// \brief R2P class for reconstructions in 3 dimensions with 2 plane (hence
/// 3D2P)
template <class CellType>
class R2P_3D2P : public R2PCommon<CellType, R2P_3D2P_columns> {
 public:
  using cell_type = CellType;
  /// \brief Number of columns for R2P_3D2P optimization.
  static constexpr UnsignedIndex_t columns_m = R2P_3D2P_columns;

  /// \brief Default constructor
  R2P_3D2P(void) = default;
  /// \brief Default destructor
  ~R2P_3D2P(void) = default;

  /// \brief Perform optimization and return planar separator.
  PlanarSeparator solve(const R2PNeighborhood<CellType>& a_neighborhood,
                        const PlanarSeparator& a_reconstruction);

  /// \brief Initialize simulation parameters that will
  /// be needed during the optimization.
  void setup(const R2PNeighborhood<CellType>& a_neighborhood,
             const PlanarSeparator& a_reconstruction);

  void allocateMatrices(const UnsignedIndex_t a_neighborhood_size);

  /// \brief Updates current guess reconstruction using a_delta
  /// and calculates the new reconstruction as well as
  /// stores weighted guess_value vector in `guess_values_m`.
  ///
  /// Delta may need to be modified as well when placing the
  /// the planes requires reprojection of the distances
  /// onto a volume-conserving distance.
  void updateGuess(Eigen::Matrix<double, columns_m, 1>* a_delta);

  /// \brief Save current guess as best and update all internally
  /// saved variables such as reference frame, beta, etc.
  void updateBestGuess(void);

  /// \brief Returns whether the minimum is reached.
  bool minimumReached(const Eigen::Matrix<double, columns_m, 1>& delta);

  /// \brief Return initial delta to be used in optimization
  auto getDefaultInitialDelta(void)
      -> const Eigen::Matrix<double, columns_m, 1>&;

  /// \brief Return best beta_m value
  const double& getBestBeta(void);

 private:
  /// \brief Return rotation for R2P dictated by elements in `a_delta`.
  ///
  /// The rotation order is:
  /// - Rotate by `a_delta(0)` radians around `a_reference_frame[0]`
  /// - Rotate by `a_delta(1)` radians around `a_reference_frame[1]`
  /// - Rotate by `a_delta(2)` radians around `a_reference_frame[2]`
  UnitQuaternion getDeltaRotationQuat(
      const ReferenceFrame& a_reference_frame,
      const Eigen::Matrix<double, columns_m, 1>& a_delta);

  /// \brief Return a reconstruction from R2P parameters
  PlanarSeparator getReconstructionFromR2PParam(
      const ReferenceFrame& a_reference_frame, const double a_beta,
      const double* a_distances);

  /// \brief Initial delta to supply to Levenberg Marquardt
  ///  to calculate Jacobian.
  Eigen::Matrix<double, columns_m, 1> initial_delta_m;
  /// \brief Current guess for beta, the rotation from shared normal
  ///  to two plane normals.
  double guess_beta_m;
  /// \brief Current guess for distance to each plane.
  double guess_distances_m[2];
  /// \brief Current best value for beta, the rotation from shared normal
  ///  to two plane normals, that gives the least error.
  double best_beta_m;
  /// \brief Current best values for distance to each plane that gives
  /// the least error.
  double best_distances_m[2];
};

/// \brief This class just calls the R2PType functions
/// but allows debug statements to be printed.
///
/// This class masks R2PType functions to allow calling
/// the function and then all printing any debug information
/// wanted. It also stores the guessed reconstructions and
/// saved "best" reconstructions in order to reprint them
/// out later. The printing format of this is to plot
/// the resulting vertices of the plane polygons
/// to be plotted by `R2P/references/r2p_history_plot.m`.
/// To do this, the information printed to screen needs to be
/// copied to `R2P/references/reconstruction_history.m`.
///
/// Requirements for R2PType class:
/// - R2PType class must meet all requirements
/// needed by the `LevenbergMarquardt` class in `optimizers.h`.
/// Right now, this is mainly R2P_2D1P, R2P_3D1P, R2P_2D2P, and R2P_3D2P.
template <class R2PType>
class R2PDebug : public R2PType {
 public:
  /// \brief Default constructor
  R2PDebug(void) = default;

  /// \brief Calls updateGuess from the base class
  /// and stores the current reconstruction.
  void updateGuess(Eigen::Matrix<double, R2PType::columns_m, 1>* a_delta);

  /// \brief Solve call for debugging R2P.
  PlanarSeparator solve(
      const R2PNeighborhood<typename R2PType::cell_type>& a_neighborhood,
      const PlanarSeparator& a_reconstruction);

  /// \brief Calls updateBestGuess from the base class
  /// and stores the current reconstruction.
  void updateBestGuess(void);

  /// \brief Shadowed getFinalReconstruction call that writes out
  /// stored best reconstructions and then returns the reconstruction.
  PlanarSeparator getFinalReconstruction(void);

  /// \brief Write the ConvexPolygons in the reconstruction out
  /// to std::cout, tagged with the given iteration number.
  void writeOutPlane(const PlanarSeparator& a_reconstruction,
                     const std::string& a_prefix,
                     const std::size_t a_iteration_number);

  /// \brief Write out the centroids and weights to
  /// to enable visualization of what optimization is driving towards.
  void writeOutCentroidsAndWeights(void);

  /// \brief Default destructor
  ~R2PDebug(void) = default;

 private:
  /// \brief Saved guess reconstructions encountered during optimization.
  std::vector<PlanarSeparator> guess_reconstruction_history;
  /// \brief Saved best reconstructions accepted during optimization.
  std::vector<PlanarSeparator> best_reconstruction_history;
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/r2p_optimization.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_OPTIMIZATION_H_
