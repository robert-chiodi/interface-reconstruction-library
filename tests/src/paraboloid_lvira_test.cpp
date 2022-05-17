// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"
#include "irl/interface_reconstruction_methods/constrained_optimization_behavior.h"
#include "irl/interface_reconstruction_methods/plvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/moments/volume_moments_with_gradient.h"
#include "irl/moments/volume_with_gradient.h"
#include "irl/optimization/constrained_levenberg_marquardt.h"
#include "irl/paraboloid_reconstruction/gradient_paraboloid.h"

#include <cmath>
#include <random>

#include "gtest/gtest.h"

#include "irl/data_structures/small_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_paraboloid.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_paraboloid.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/parametrized_surface.h"
#include "irl/planar_reconstruction/planar_separator.h"

#include <Eigen/Dense>  // Eigen header
#include "tests/src/basic_mesh.h"
#include "tests/src/data.h"
#include "tests/src/vtk.h"

namespace IRL {

TEST(ParaboloidLVIRA, ProgressiveDistanceSolver) {
  AlignedParaboloid aligned_paraboloid;
  aligned_paraboloid.a() = 10.0;
  aligned_paraboloid.b() = 20.0;
  std::array<double, 3> angles{{M_PI / 10.0, M_PI / 5.0}};
  // std::array<double, 3> angles{{0.0, 0.0}};
  std::array<double, 3> translations{{0.0, 0.0, 0.0}};
  ReferenceFrame frame(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                       Normal(0.0, 0.0, 1.0));
  UnitQuaternion x_rotation(angles[0], frame[0]);
  UnitQuaternion y_rotation(angles[1], frame[1]);
  UnitQuaternion z_rotation(angles[2], frame[2]);
  auto total_rotation = x_rotation * y_rotation * z_rotation;
  total_rotation.normalize();
  frame = total_rotation * frame;
  auto datum = -Pt::fromArray(translations);
  Paraboloid paraboloid(datum, frame, aligned_paraboloid.a(),
                        aligned_paraboloid.b());

  auto cube = RectangularCuboid::fromBoundingPts(Pt(-1.0, -1.0, -1.0),
                                                 Pt(1.0, 1.0, 1.0));
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  std::uniform_real_distribution<double> random_vfrac(0.0, 1.0);

  double max_error = 0.0;
  double tolerance = 1.0e-14;
  int Ntest = 1e2;
  for (int i = 0; i < Ntest; ++i) {
    double vfrac_required = random_vfrac(eng);
    ProgressiveDistanceSolverParaboloid<RectangularCuboid> solver(
        cube, vfrac_required, tolerance, paraboloid);
    double distance = solver.getDistance();
    Paraboloid new_paraboloid(datum + distance * frame[2], frame,
                              aligned_paraboloid.a(), aligned_paraboloid.b());
    double error =
        std::fabs(getVolumeFraction(cube, new_paraboloid) - vfrac_required);
    max_error = std::max(max_error, error);
    // if (error > 1.0e-12) {
    // std::cout << std::setprecision(20) << "VFRAC = " << vfrac_required;
    // std::cout << std::setprecision(20) << " -- distance = " << distance;
    // std::cout << std::setprecision(20) << " -- error = " << error <<
    // std::endl;
    if (error > 1.0e-12) {
      exit(-1);
    }
  }

  std::cout << std::setprecision(20) << Ntest
            << " tests -- max error = " << max_error << std::endl;

  EXPECT_NEAR(max_error, 0.0, 10.0 * tolerance);
}

// Simple class to look for zeros on with optimization
template <class CellType, UnsignedIndex_t kMeasures,
          UnsignedIndex_t kParameters, UnsignedIndex_t kConstraints>
class PLVIRA_test {
  using cell_type = CellType;

 public:
  PLVIRA_test() = default;

  PLVIRA_test(const PLVIRANeighborhood<CellType>* a_neighborhood,
              const Paraboloid a_reconstruction) {
    neighborhood_m = a_neighborhood;
    ref_volume_m = neighborhood_m->getConstrainedCell().calculateVolume();
    ref_length_m = std::pow(ref_volume_m, 1.0 / 3.0);
    ref_moment_m = ref_length_m * ref_volume_m;
    best_reconstruction_m = a_reconstruction;
    for (UnsignedIndex_t i = 0; i < kMeasures; ++i) {
      correct_values_m(i) = neighborhood_m->getStoredMoments(i) / ref_volume_m;
      weights_m.diagonal()[i] = neighborhood_m->getWeight(i);
    }
    for (UnsignedIndex_t i = kMeasures; i < kMeasures + kConstraints; ++i) {
      weights_m.diagonal()[i] = 1.0;
    }
    correct_constraints_m(0) =
        neighborhood_m->getConstraints().volume() / ref_volume_m;
    correct_constraints_m(1) =
        neighborhood_m->getConstraints().centroid()[0] / ref_moment_m;
    correct_constraints_m(2) =
        neighborhood_m->getConstraints().centroid()[1] / ref_moment_m;
    correct_constraints_m(3) =
        neighborhood_m->getConstraints().centroid()[2] / ref_moment_m;
    // correct_constraints_m(0) = 0.0;
    // correct_constraints_m(1) = 0.0;
    // correct_constraints_m(2) = 0.0;
    // correct_constraints_m(3) = 0.0;
    for (UnsignedIndex_t i = 0; i < kConstraints; ++i) {
      correct_values_m(kMeasures + i) = correct_constraints_m(i);
    }
  }

  Eigen::DiagonalMatrix<double, kMeasures + kConstraints>& getWeights(void) {
    return weights_m;
  }

  bool errorTooHigh(const double a_error) {
    return a_error >
           optimization_behavior_m.acceptable_squared_constraint_error;
  }

  bool subErrorTooHigh(const double a_error) {
    return a_error > optimization_behavior_m.acceptable_squared_error;
  }

  bool maxErrorTooHigh(Eigen::Matrix<double, kParameters, 1>* a_error) {
    double linf_rhs = 0.0;
    for (UnsignedIndex_t i = 0; i < kParameters; ++i) {
      linf_rhs = std::max(linf_rhs, std::fabs((*a_error)(i)));
    }
    return linf_rhs > optimization_behavior_m.acceptable_max_error;
  }

  bool iterationTooHigh(const UnsignedIndex_t a_iterations) {
    return a_iterations > optimization_behavior_m.maximum_iterations;
  }

  bool subIterationTooHigh(const UnsignedIndex_t a_sub_iterations) {
    return a_sub_iterations > optimization_behavior_m.maximum_sub_iterations;
  }

  double getDampingParameterInitialFactor(void) {
    return optimization_behavior_m.damping_param_initial_factor;
  }

  double getDampingParameterIncreaseFactor(void) {
    return optimization_behavior_m.damping_param_increase_factor;
  }

  double getPenaltyParameterIncreaseFactor(void) {
    return optimization_behavior_m.penalty_param_increase_factor;
  }

  void updateGuess(const Eigen::Matrix<double, kParameters, 1>* a_delta,
                   const double a_penalty,
                   const Eigen::Matrix<double, kConstraints, 1> a_multipliers) {
    // Get old datum, frame, and curvatures
    const auto& datum = best_reconstruction_m.getDatum();
    const auto& frame = best_reconstruction_m.getReferenceFrame();
    const auto& algnd_para = best_reconstruction_m.getAlignedParaboloid();
    // Create quaternion for rotating the frame
    UnitQuaternion x_rotation((*a_delta)(3), frame[0]);
    UnitQuaternion y_rotation((*a_delta)(4), frame[1]);
    UnitQuaternion z_rotation((*a_delta)(5), frame[2]);
    auto total_rotation = x_rotation * y_rotation * z_rotation;
    total_rotation.normalize();
    // Shift datum
    const Pt new_datum = datum + (*a_delta)(0) * frame[0] +
                         (*a_delta)(1) * frame[1] + (*a_delta)(2) * frame[2];
    // Rotate frame
    const ReferenceFrame new_frame = total_rotation * frame;
    // Update curvatures
    const double new_a = algnd_para.a() + (*a_delta)(6);
    const double new_b = algnd_para.b() + (*a_delta)(7);
    // Return new paraboloid guess
    guess_reconstruction_m = Paraboloid(new_datum, new_frame, new_a, new_b);
    guess_values_m.setZero();
    for (UnsignedIndex_t i = 0; i < kMeasures; ++i) {
      guess_values_m(i) = getVolumeMoments<Volume>(neighborhood_m->getCell(i),
                                                   guess_reconstruction_m) /
                          ref_volume_m;
    }
    auto moments_constrained_cell = getVolumeMoments<VolumeMoments>(
        neighborhood_m->getConstrainedCell(), guess_reconstruction_m);
    guess_values_m(kMeasures) =
        moments_constrained_cell.volume() / ref_volume_m;
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      guess_values_m(kMeasures + 1 + d) =
          moments_constrained_cell.centroid()[d] / ref_moment_m;
    }
  }

  void updateBestGuess(void) {
    best_values_m = guess_values_m;
    best_reconstruction_m = guess_reconstruction_m;
  }

  void storeBestGuess(void) {
    stored_best_reconstruction_m = best_reconstruction_m;
  }

  void restoreBestGuess(void) {
    best_reconstruction_m = stored_best_reconstruction_m;
  }

  void storeInitialGuess(void) {
    stored_initial_reconstruction_m = best_reconstruction_m;
  }

  void restoreInitialGuess(void) {
    best_reconstruction_m = stored_initial_reconstruction_m;
  }

  double calculateScalarError(
      const double a_penalty,
      const Eigen::Matrix<double, kConstraints, 1> a_multipliers) {
    double error = 0.0;
    for (UnsignedIndex_t i = 0; i < kMeasures; ++i) {
      error += (guess_values_m(i) - correct_values_m(i)) *
               (guess_values_m(i) - correct_values_m(i)) *
               weights_m.diagonal()[i];
    }
    for (UnsignedIndex_t i = 0; i < kConstraints; ++i) {
      error +=
          a_penalty *
          ((guess_values_m(kMeasures + i) - correct_values_m(kMeasures + i)) +
           a_multipliers(i) / (2.0 * a_penalty)) *
          ((guess_values_m(kMeasures + i) - correct_values_m(kMeasures + i)) +
           a_multipliers(i) / (2.0 * a_penalty)) *
          weights_m.diagonal()[kMeasures + i];
    }
    return error;
  }

  void clipChange(Eigen::Matrix<double, kParameters, 1>* a_delta) {
    double scale = 1.0;
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      if (std::fabs((*a_delta)(d)) > 0.5 * ref_length_m) {
        scale = std::min(scale, 0.5 * ref_length_m / std::fabs((*a_delta)(d)));
      }
      if (std::fabs((*a_delta)(3 + d)) > 0.1 * M_PI) {
        scale = std::min(scale, 0.1 * M_PI / std::fabs((*a_delta)(3 + d)));
      }
    }
    if (std::fabs((*a_delta)(6)) > 0.1 / ref_length_m) {
      scale = std::min(scale, 0.1 / ref_length_m / std::fabs((*a_delta)(6)));
    }
    if (std::fabs((*a_delta)(7)) > 0.1 / ref_length_m) {
      scale = std::min(scale, 0.1 / ref_length_m / std::fabs((*a_delta)(7)));
    }
    if (scale != 1.0) {
      for (UnsignedIndex_t d = 0; d < kParameters; ++d) {
        (*a_delta)(d) *= scale;
      }
    }
  }

  // void calculateConstraintVectorError(
  //     Eigen::Matrix<double, kConstraints, 1>* a_constraints) {
  //   assert(neighborhood_m->size() == kMeasures);
  //   (*a_constraints).setZero();
  //   const auto& constrained_cell = neighborhood_m->getConstrainedCell();
  //   const auto constrained_moments = getVolumeMoments<VolumeMoments>(
  //       constrained_cell, best_reconstruction_m);
  //   (*a_constraints)(0) =
  //       (correct_constraints_m(0) - constrained_moments.volume()) /
  //       ref_volume_m;
  //   (*a_constraints)(1) = =
  //       correct_constraints_m(1, 0) - moments_with_gradients.centroid()[0];
  //   (*a_constraints)(2) =
  //       correct_constraints_m(2, 0) - moments_with_gradients.centroid()[1];
  //   (*a_constraints)(3) =
  //       correct_constraints_m(3, 0) - moments_with_gradients.centroid()[2];
  // }

  void calculateVectorErrorAndJacobian(
      Eigen::Matrix<double, kMeasures + kConstraints, 1>* a_error_vector,
      Eigen::Matrix<double, kParameters, kMeasures + kConstraints>*
          a_jacobian_transpose,
      Eigen::Matrix<double, kConstraints, 1>* a_constraints,
      const double a_penalty,
      Eigen::Matrix<double, kConstraints, 1> a_multipliers) {
    assert(neighborhood_m->size() == kMeasures);
    (*a_error_vector).setZero();
    (*a_constraints).setZero();
    (*a_jacobian_transpose).setZero();
    using MyGradientType = ParaboloidGradientLocal;
    using MyPtType = PtWithGradient<MyGradientType>;
    // Jacobian of unconstrained measure
    for (UnsignedIndex_t i = 0; i < neighborhood_m->size(); ++i) {
      const auto& cell = neighborhood_m->getCell(i);
      const auto cell_with_data =
          StoredRectangularCuboid<MyPtType>::fromOtherPolytope(cell);
      const auto moments_with_gradients =
          getVolumeMoments<VolumeWithGradient<MyGradientType>>(
              cell_with_data, best_reconstruction_m);
      const auto& gradient = moments_with_gradients.volume_gradient();
      const double vf = moments_with_gradients.volume() / ref_volume_m;
      if ((vf < global_constants::VF_LOW &&
           correct_values_m(i) < global_constants::VF_LOW) ||
          (vf > global_constants::VF_HIGH &&
           correct_values_m(i) > global_constants::VF_HIGH)) {
        // We are good
      } else {
        (*a_error_vector)(i) = (correct_values_m(i) -
                                moments_with_gradients.volume() / ref_volume_m);
        (*a_jacobian_transpose)(0, i) = gradient.getGradTx() / ref_volume_m;
        (*a_jacobian_transpose)(1, i) = gradient.getGradTy() / ref_volume_m;
        (*a_jacobian_transpose)(2, i) = gradient.getGradTz() / ref_volume_m;
        (*a_jacobian_transpose)(3, i) = gradient.getGradRx() / ref_volume_m;
        (*a_jacobian_transpose)(4, i) = gradient.getGradRy() / ref_volume_m;
        (*a_jacobian_transpose)(5, i) = gradient.getGradRz() / ref_volume_m;
        (*a_jacobian_transpose)(6, i) = gradient.getGradA() / ref_volume_m;
        (*a_jacobian_transpose)(7, i) = gradient.getGradB() / ref_volume_m;
      }
    }
    // Jacobian of constraints
    const auto& constrained_cell = neighborhood_m->getConstrainedCell();
    const auto constrained_moments = getVolumeMoments<VolumeMoments>(
        constrained_cell, best_reconstruction_m);
    const auto constrained_cell_with_data =
        StoredRectangularCuboid<MyPtType>::fromOtherPolytope(constrained_cell);
    const auto constrained_moments_with_gradients =
        getVolumeMoments<VolumeMomentsWithGradient<MyGradientType>>(
            constrained_cell_with_data, best_reconstruction_m);

    // if (constrained_moments_with_gradients.volume() / ref_volume_m < -100.0)
    // {
    //   exit(-1);
    // }

    const auto& constrained_volume_gradient =
        constrained_moments_with_gradients.volume_gradient();
    const auto& constrained_centroid_gradient =
        constrained_moments_with_gradients.centroid().getData();
    (*a_constraints)(0) =
        correct_constraints_m(0) -
        constrained_moments_with_gradients.volume() / ref_volume_m;
    (*a_jacobian_transpose)(0, kMeasures) =
        constrained_volume_gradient.getGradTx() / ref_volume_m;
    (*a_jacobian_transpose)(1, kMeasures) =
        constrained_volume_gradient.getGradTy() / ref_volume_m;
    (*a_jacobian_transpose)(2, kMeasures) =
        constrained_volume_gradient.getGradTz() / ref_volume_m;
    (*a_jacobian_transpose)(3, kMeasures) =
        constrained_volume_gradient.getGradRx() / ref_volume_m;
    (*a_jacobian_transpose)(4, kMeasures) =
        constrained_volume_gradient.getGradRy() / ref_volume_m;
    (*a_jacobian_transpose)(5, kMeasures) =
        constrained_volume_gradient.getGradRz() / ref_volume_m;
    (*a_jacobian_transpose)(6, kMeasures) =
        constrained_volume_gradient.getGradA() / ref_volume_m;
    (*a_jacobian_transpose)(7, kMeasures) =
        constrained_volume_gradient.getGradB() / ref_volume_m;
    // double norm_grad = 0.0;
    // for (int i = 0; i < 8; i++) {
    //   norm_grad += (constrained_volume_gradient.getGrad()[i] +
    //                 constrained_centroid_gradient[0].getGrad()[i] +
    //                 constrained_centroid_gradient[1].getGrad()[i] +
    //                 constrained_centroid_gradient[2].getGrad()[i]) *
    //                (constrained_volume_gradient.getGrad()[i] +
    //                 constrained_centroid_gradient[0].getGrad()[i] +
    //                 constrained_centroid_gradient[1].getGrad()[i] +
    //                 constrained_centroid_gradient[2].getGrad()[i]);
    // }
    // std::cout << "Constrained grad norm = " << std::sqrt(norm_grad)
    //           << std::endl;
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      (*a_constraints)(1 + d) =
          correct_constraints_m(1 + d, 0) -
          constrained_moments_with_gradients.centroid().getPt()[d] /
              ref_moment_m;
      (*a_jacobian_transpose)(0, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradTx() / ref_moment_m;
      (*a_jacobian_transpose)(1, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradTy() / ref_moment_m;
      (*a_jacobian_transpose)(2, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradTz() / ref_moment_m;
      (*a_jacobian_transpose)(3, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradRx() / ref_moment_m;
      (*a_jacobian_transpose)(4, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradRy() / ref_moment_m;
      (*a_jacobian_transpose)(5, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradRz() / ref_moment_m;
      (*a_jacobian_transpose)(6, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradA() / ref_moment_m;
      (*a_jacobian_transpose)(7, kMeasures + 1 + d) =
          constrained_centroid_gradient[d].getGradB() / ref_moment_m;
    }
    for (UnsignedIndex_t i = 0; i < kConstraints; ++i) {
      (*a_error_vector)(kMeasures + i) =
          std::sqrt(a_penalty) *
          ((*a_constraints)(i) + a_multipliers(i) / (2.0 * a_penalty));
      for (UnsignedIndex_t j = 0; j < kParameters; ++j) {
        (*a_jacobian_transpose)(j, kMeasures + i) *= std::sqrt(a_penalty);
      }
    }
  }

  bool minimumReached(const Eigen::Matrix<double, kParameters, 1> a_delta) {
    return a_delta.squaredNorm() < 1.0e-16;
  }

  Paraboloid& getBestReconstruction(void) { return best_reconstruction_m; }

  ~PLVIRA_test() = default;

 private:
  const PLVIRANeighborhood<CellType>* neighborhood_m;
  Eigen::DiagonalMatrix<double, kMeasures + kConstraints> weights_m;
  Eigen::Matrix<double, kMeasures + kConstraints, 1> correct_values_m;
  Eigen::Matrix<double, kConstraints, 1> correct_constraints_m;
  double ref_length_m;
  double ref_volume_m;
  double ref_moment_m;
  ConstrainedOptimizationBehavior optimization_behavior_m;
  Eigen::Matrix<double, kMeasures + kConstraints, 1> guess_values_m;
  Eigen::Matrix<double, kMeasures + kConstraints, 1> best_values_m;
  Paraboloid guess_reconstruction_m;
  Paraboloid best_reconstruction_m;
  Paraboloid stored_best_reconstruction_m;
  Paraboloid stored_initial_reconstruction_m;
};  // namespace IRL

TEST(ParaboloidLVIRA, PLVIRA) {
  using PLVIRA_type = PLVIRA_test<RectangularCuboid, 26, 8, 4>;

  // Construct NxNxN mesh
  int n_cells = 10;
  constexpr const int number_of_ghost_cells = 0;
  const double dx = 1.0 / static_cast<double>(n_cells);
  Pt lower_domain(-0.5, -0.5, -0.5);
  Pt upper_domain(0.5, 0.5, 0.5);
  BasicMesh mesh(n_cells, n_cells, n_cells, number_of_ghost_cells);
  mesh.setCellBoundaries(lower_domain, upper_domain);

  // Initialize VF from sphere
  auto center = Pt(0.001, 0.001, 0.001);
  double radius = 0.2;
  double curvature = 1.0 / radius;

  Data<double> liquid_vf(mesh);
  Data<VolumeMoments> liquid_moments(mesh);
  Data<IRL::Paraboloid> liquid_gas_interface(mesh);

  int sub_div = 10;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        double sub_dx =
            (mesh.x(i + 1) - mesh.x(i)) / static_cast<double>(sub_div);
        double sub_dy =
            (mesh.y(j + 1) - mesh.y(j)) / static_cast<double>(sub_div);
        double sub_dz =
            (mesh.z(k + 1) - mesh.z(k)) / static_cast<double>(sub_div);
        liquid_moments(i, j, k) = VolumeMoments::fromScalarConstant(0.0);
        for (int kk = 0; kk < sub_div; ++kk) {
          for (int jj = 0; jj < sub_div; ++jj) {
            for (int ii = 0; ii < sub_div; ++ii) {
              const Pt lower_cell_pt(
                  mesh.x(i) + static_cast<double>(ii) * sub_dx,
                  mesh.y(j) + static_cast<double>(jj) * sub_dy,
                  mesh.z(k) + static_cast<double>(kk) * sub_dz);
              const Pt upper_cell_pt(
                  mesh.x(i) + static_cast<double>(ii + 1) * sub_dx,
                  mesh.y(j) + static_cast<double>(jj + 1) * sub_dy,
                  mesh.z(k) + static_cast<double>(kk + 1) * sub_dz);
              const auto sub_cell = RectangularCuboid::fromBoundingPts(
                  lower_cell_pt, upper_cell_pt);
              Normal normal = sub_cell.calculateCentroid() - center;
              normal.normalize();
              int largest_dir = 0;
              if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
                largest_dir = 1;
              if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
                largest_dir = 2;
              ReferenceFrame frame;
              if (largest_dir == 0)
                frame[0] = crossProduct(normal, Normal(0.0, 1.0, 0.0));
              else if (largest_dir == 0)
                frame[0] = crossProduct(normal, Normal(0.0, 0.0, 1.0));
              else
                frame[0] = crossProduct(normal, Normal(1.0, 0.0, 0.0));
              frame[0].normalize();
              frame[1] = crossProduct(normal, frame[0]);
              frame[2] = normal;
              Paraboloid paraboloid(Pt(center + radius * normal), frame,
                                    0.5 * curvature, 0.5 * curvature);
              liquid_moments(i, j, k) +=
                  getVolumeMoments<VolumeMoments>(sub_cell, paraboloid);
            }
          }
        }
        {
          const Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          liquid_vf(i, j, k) =
              liquid_moments(i, j, k).volume() / cell.calculateVolume();
          liquid_vf(i, j, k) = std::max(liquid_vf(i, j, k), 0.0);
          liquid_vf(i, j, k) = std::min(liquid_vf(i, j, k), 1.0);

          Normal normal = cell.calculateCentroid() - center;
          normal.normalize();
          if (i == 1 && j == 1 & k == 1) normal[0] += 0.0;
          normal.normalize();
          int largest_dir = 0;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
            largest_dir = 1;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
            largest_dir = 2;
          ReferenceFrame frame;
          if (largest_dir == 0)
            frame[0] = crossProduct(normal, Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 0)
            frame[0] = crossProduct(normal, Normal(0.0, 0.0, 1.0));
          else
            frame[0] = crossProduct(normal, Normal(1.0, 0.0, 0.0));
          frame[0].normalize();
          frame[1] = crossProduct(normal, frame[0]);
          frame[2] = normal;
          double a = 0.5 * curvature;
          double b = 0.5 * curvature;
          if (i == 1 && j == 1 & k == 1) {
            a = 0.5 * curvature;
            b = 0.5 * curvature;
          }
          // a = 0.0;
          // b = 0.0;
          Paraboloid paraboloid(center + radius * normal, frame, a, b);
          ProgressiveDistanceSolverParaboloid<RectangularCuboid> solver(
              cell, liquid_vf(i, j, k), 1.0e-14, paraboloid);
          Paraboloid new_paraboloid(
              Pt(center + (radius + solver.getDistance()) * normal), frame, a,
              b);
          liquid_gas_interface(i, j, k) = new_paraboloid;
        }
      }
    }
  }

  // Initialize folders/mesh for very simple I/O
  int viz_output = 0;
  double time = 0.0;
  VTKOutput vtk_io("viz_out", "viz", mesh);
  vtk_io.addData("VOF", liquid_vf);
  vtk_io.writeVTKFile(time);
  std::vector<ParametrizedSurfaceOutput> surfaces;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (liquid_vf(i, j, k) >= global_constants::VF_LOW &&
            liquid_vf(i, j, k) <= global_constants::VF_HIGH) {
          const Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          auto volume_and_surface = getVolumeMoments<
              AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
              cell, liquid_gas_interface(i, j, k));
          auto& surface = volume_and_surface.getSurface();
          surface.setLengthScale(std::pow(cell.calculateVolume(), 1.0 / 3.0) /
                                 10.0);
          surfaces.push_back(surface);
        }
      }
    }
  }
  vtk_io.writeVTKInterface(time, surfaces, true);
  ++viz_output;

  surfaces.clear();
  for (int i = mesh.imin() + 1; i <= mesh.imax() - 1; ++i) {
    for (int j = mesh.jmin() + 1; j <= mesh.jmax() - 1; ++j) {
      for (int k = mesh.kmin() + 1; k <= mesh.kmax() - 1; ++k) {
        if (liquid_vf(i, j, k) >= global_constants::VF_LOW &&
            liquid_vf(i, j, k) <= global_constants::VF_HIGH) {
          std::cout << "Creating PLVIRA neighborhood " << std::endl;
          // Intialize LVIRA 3x3x3 neighbourhood
          Paraboloid paraboloid_guess = liquid_gas_interface(i, j, k);
          RectangularCuboid stencil_cells[27];
          std::array<double, 27> stencil_weights;
          double cellVF[27];
          PLVIRANeighborhood<RectangularCuboid> neighborhood_VF;
          for (int kk = 0; kk < 3; ++kk) {
            for (int jj = 0; jj < 3; ++jj) {
              for (int ii = 0; ii < 3; ++ii) {
                const Pt lower_cell_pt(mesh.x(i + ii - 1), mesh.y(j + jj - 1),
                                       mesh.z(k + kk - 1));
                const Pt upper_cell_pt(mesh.x(i + ii), mesh.y(j + jj),
                                       mesh.z(k + kk));
                stencil_cells[ii + jj * 3 + kk * 9] =
                    RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                       upper_cell_pt);
                cellVF[ii + jj * 3 + kk * 9] =
                    liquid_moments(i + ii - 1, j + jj - 1, k + kk - 1).volume();

                double weight = 1.0 / static_cast<double>(1 + std::abs(ii - 1) +
                                                          std::abs(jj - 1) +
                                                          std::abs(kk - 1));
                // weight *= 2.0 * std::fabs(liquid_vf(i + ii - 1, j + jj - 1,
                //                                     k + kk - 1) -
                //                           0.5);
                if (ii == 1 && jj == 1 && kk == 1) {
                  neighborhood_VF.setConstrainedCell(
                      &stencil_cells[ii + jj * 3 + kk * 9],
                      liquid_moments(i + ii - 1, j + jj - 1, k + kk - 1));
                } else {
                  neighborhood_VF.addMember(
                      &stencil_cells[ii + jj * 3 + kk * 9],
                      &cellVF[ii + jj * 3 + kk * 9], weight);
                }
              }
            }
          }

          PLVIRA_type plvira_object(&neighborhood_VF, paraboloid_guess);
          ConstrainedLevenbergMarquardt<PLVIRA_type, 26, 8, 4> solver;
          std::cout << "Solving... " << std::endl;
          solver.solve(&plvira_object);
          std::cout << "Solved in   : " << solver.getIterationCount()
                    << " multiplier updates" << std::endl;
          std::cout << "              " << solver.getSubIterationCount()
                    << " sub-iterations" << std::endl;
          std::cout << " with error : " << std::endl
                    << solver.getConstraintErrorVector() << std::endl;

          // PLVIRA_type plvira_object_refinement(
          //     &neighborhood_VF, plvira_object.getBestReconstruction());
          // ConstrainedLevenbergMarquardt<PLVIRA_type, 26, 8, 4>
          //     solver_refinement;
          // std::cout << "Refining... " << std::endl;
          // solver_refinement.solve(&plvira_object);
          // std::cout << "Refined in   : "
          //           << solver_refinement.getIterationCount()
          //           << " multiplier updates" << std::endl;
          // std::cout << "              "
          //           << solver_refinement.getSubIterationCount()
          //           << " sub-iterations" << std::endl;
          // std::cout << " with error : " << std::endl
          //           << solver_refinement.getConstraintErrorVector()
          //           << std::endl;

          const auto cell = neighborhood_VF.getConstrainedCell();
          auto volume_and_surface = getVolumeMoments<
              AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
              cell, plvira_object.getBestReconstruction());
          auto& surface = volume_and_surface.getSurface();
          surface.setLengthScale(std::pow(cell.calculateVolume(), 1.0 / 3.0) /
                                 30.0);
          surfaces.push_back(surface);
        }
      }
    }
  }
  std::cout << "Printing final solution " << std::endl;
  vtk_io.writeVTKInterface(time, surfaces, true);
}

void ComputeJacobianAndF(Paraboloid& paraboloid,
                         Eigen::Matrix<double, 8, 27>& jacobian_transpose,
                         Eigen::Matrix<double, 27, 1>& f,
                         LVIRANeighborhood<RectangularCuboid> lvira_stencil);

TEST(ParaboloidLVIRA, LVIRA) {
  // Construct NxNxN mesh
  int n_cells = 3;
  constexpr const int number_of_ghost_cells = 0;
  const double dx = 1.0 / static_cast<double>(n_cells);
  Pt lower_domain(-0.5, -0.5, -0.5);
  Pt upper_domain(0.5, 0.5, 0.5);
  BasicMesh mesh(n_cells, n_cells, n_cells, number_of_ghost_cells);
  mesh.setCellBoundaries(lower_domain, upper_domain);

  // Initialize VF from sphere
  auto center = Pt(-0.4, -0.5, -0.52);
  double radius = 0.5 * std::sqrt(3.0);
  double curvature = 1.0 / radius;

  Data<double> liquid_volume_fraction(mesh);
  Data<IRL::Paraboloid> liquid_gas_interface(mesh);

  int sub_div = 10;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        double sub_dx =
            (mesh.x(i + 1) - mesh.x(i)) / static_cast<double>(sub_div);
        double sub_dy =
            (mesh.y(j + 1) - mesh.y(j)) / static_cast<double>(sub_div);
        double sub_dz =
            (mesh.z(k + 1) - mesh.z(k)) / static_cast<double>(sub_div);
        liquid_volume_fraction(i, j, k) = 0.0;
        for (int kk = 0; kk < sub_div; ++kk) {
          for (int jj = 0; jj < sub_div; ++jj) {
            for (int ii = 0; ii < sub_div; ++ii) {
              const Pt lower_cell_pt(
                  mesh.x(i) + static_cast<double>(ii) * sub_dx,
                  mesh.y(j) + static_cast<double>(jj) * sub_dy,
                  mesh.z(k) + static_cast<double>(kk) * sub_dz);
              const Pt upper_cell_pt(
                  mesh.x(i) + static_cast<double>(ii + 1) * sub_dx,
                  mesh.y(j) + static_cast<double>(jj + 1) * sub_dy,
                  mesh.z(k) + static_cast<double>(kk + 1) * sub_dz);
              const auto sub_cell = RectangularCuboid::fromBoundingPts(
                  lower_cell_pt, upper_cell_pt);
              Normal normal = sub_cell.calculateCentroid() - center;
              normal.normalize();
              int largest_dir = 0;
              if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
                largest_dir = 1;
              if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
                largest_dir = 2;
              ReferenceFrame frame;
              if (largest_dir == 0)
                frame[0] = crossProduct(normal, Normal(0.0, 1.0, 0.0));
              else if (largest_dir == 0)
                frame[0] = crossProduct(normal, Normal(0.0, 0.0, 1.0));
              else
                frame[0] = crossProduct(normal, Normal(1.0, 0.0, 0.0));
              frame[0].normalize();
              frame[1] = crossProduct(normal, frame[0]);
              frame[2] = normal;
              Paraboloid paraboloid(Pt(center + radius * normal), frame,
                                    0.5 * curvature, 0.5 * curvature);
              liquid_volume_fraction(i, j, k) +=
                  getVolumeFraction(sub_cell, paraboloid);
            }
          }
        }
        {
          liquid_volume_fraction(i, j, k) /=
              static_cast<double>(sub_div * sub_div * sub_div);
          liquid_volume_fraction(i, j, k) =
              std::min(liquid_volume_fraction(i, j, k), 1.0);
          liquid_volume_fraction(i, j, k) =
              std::max(liquid_volume_fraction(i, j, k), 0.0);

          const Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          Normal normal = cell.calculateCentroid() - center;
          normal.normalize();
          if (i == 1 && j == 1 & k == 1) normal[0] += 10.0;
          normal.normalize();
          int largest_dir = 0;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
            largest_dir = 1;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
            largest_dir = 2;
          ReferenceFrame frame;
          if (largest_dir == 0)
            frame[0] = crossProduct(normal, Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 0)
            frame[0] = crossProduct(normal, Normal(0.0, 0.0, 1.0));
          else
            frame[0] = crossProduct(normal, Normal(1.0, 0.0, 0.0));
          frame[0].normalize();
          frame[1] = crossProduct(normal, frame[0]);
          frame[2] = normal;
          Paraboloid paraboloid(center + radius * normal, frame,
                                0.5 * curvature, 0.5 * curvature);
          ProgressiveDistanceSolverParaboloid<RectangularCuboid> solver(
              cell, liquid_volume_fraction(i, j, k), 1.0e-14, paraboloid);
          Paraboloid new_paraboloid(
              Pt(center + (radius + solver.getDistance()) * normal), frame,
              0.5 * curvature, 0.5 * curvature);
          liquid_gas_interface(i, j, k) = new_paraboloid;
        }
      }
    }
  }

  // Initialize folders/mesh for very simple I/O
  int viz_output = 0;
  double time = 0.0;
  VTKOutput vtk_io("viz_out", "viz", mesh);
  vtk_io.addData("VOF", liquid_volume_fraction);
  vtk_io.writeVTKFile(time);
  std::vector<ParametrizedSurfaceOutput> surfaces;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (liquid_volume_fraction(i, j, k) >= global_constants::VF_LOW &&
            liquid_volume_fraction(i, j, k) <= global_constants::VF_HIGH) {
          const Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          auto volume_and_surface = getVolumeMoments<
              AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
              cell, liquid_gas_interface(i, j, k));
          auto& surface = volume_and_surface.getSurface();
          surface.setLengthScale(std::pow(cell.calculateVolume(), 1.0 / 3.0) /
                                 10.0);
          surfaces.push_back(surface);
        }
      }
    }
  }
  vtk_io.writeVTKInterface(time, surfaces, true);
  ++viz_output;

  // Intialize LVIRA 3x3x3 neighbourhood
  Paraboloid paraboloid_to_fit = liquid_gas_interface(1, 1, 1);
  RectangularCuboid stencil_cells[27];
  std::array<double, 27> stencil_weights;
  double cellVF[27];
  LVIRANeighborhood<RectangularCuboid> neighborhood_VF;
  neighborhood_VF.resize(27);
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        const Pt lower_cell_pt(mesh.x(n_cells / 2 + i - 1),
                               mesh.y(n_cells / 2 + j - 1),
                               mesh.z(n_cells / 2 + k - 1));
        const Pt upper_cell_pt(mesh.x(n_cells / 2 + i), mesh.y(n_cells / 2 + j),
                               mesh.z(n_cells / 2 + k));
        stencil_cells[i + j * 3 + k * 9] =
            RectangularCuboid::fromBoundingPts(lower_cell_pt, upper_cell_pt);
        cellVF[i + j * 3 + k * 9] = liquid_volume_fraction(
            n_cells / 2 + i - 1, n_cells / 2 + j - 1, n_cells / 2 + k - 1);
        neighborhood_VF.setMember(i + j * 3 + k * 9,
                                  &stencil_cells[i + j * 3 + k * 9],
                                  &cellVF[i + j * 3 + k * 9]);
        // stencil_weights[i + j * 3 + k * 9] = 1.0;
        // if (i == 1 && j == 1 && k == 1) {
        //   stencil_weights[i + j * 3 + k * 9] = 1.0;
        // }
        // if ((i - 1) * (j - 1) * (k - 1) == 0) {
        //   stencil_weights[i + j * 3 + k * 9] = 1.0;
        // }
        stencil_weights[i + j * 3 + k * 9] =
            1.0 / static_cast<double>(1 + std::abs(i - 1) * std::abs(i - 1) +
                                      std::abs(j - 1) * std::abs(j - 1) +
                                      std::abs(k - 1) * std::abs(k - 1));
      }
    }
  }
  neighborhood_VF.setCenterOfStencil(13);

  // for (int n = 0; n < 10; n++) {
  // Update normal
  //      {
  //        auto central_cell = neighborhood_VF.getCenterCell();
  //        auto central_moment =
  //        neighborhood_VF.getCenterCellStoredMoments(); auto
  //        volume_and_surface =
  //            getVolumeMoments<AddSurfaceOutput<Volume,
  //    ParametrizedSurfaceOutput>>( central_cell, paraboloid_to_fit); auto&
  //    surface = volume_and_surface.getSurface(); auto new_normal_local =
  //    surface.getAverageNormal(); auto& ref_frame =
  //    paraboloid_to_fit.getReferenceFrame(); auto new_normal = Normal(0.0,
  //    0.0, 0.0); for (UnsignedIndex_t d = 0; d < 3; ++d) { for
  //    (UnsignedIndex_t n = 0; n < 3; ++n) { new_normal[n] += ref_frame[d][n]
  //    * new_normal_local[d];
  //          }
  //        }
  //        int largest_dir = 0;
  //        if (std::fabs(new_normal[largest_dir]) < std::fabs(new_normal[1]))
  //          largest_dir = 1;
  //        if (std::fabs(new_normal[largest_dir]) < std::fabs(new_normal[2]))
  //          largest_dir = 2;
  //        ReferenceFrame frame;
  //        if (largest_dir == 0)
  //          frame[0] = crossProduct(new_normal, Normal(0.0, 1.0, 0.0));
  //        else if (largest_dir == 0)
  //          frame[0] = crossProduct(new_normal, Normal(0.0, 0.0, 1.0));
  //        else
  //          frame[0] = crossProduct(new_normal, Normal(1.0, 0.0, 0.0));
  //        frame[0].normalize();
  //        frame[1] = crossProduct(new_normal, frame[0]);
  //        frame[2] = new_normal;
  //        auto& aligned_paraboloid =
  //        paraboloid_to_fit.getAlignedParaboloid(); auto& datum =
  //        paraboloid_to_fit.getDatum(); auto new_paraboloid =
  //        Paraboloid(datum, frame, aligned_paraboloid.a(),
  //                                         aligned_paraboloid.b());
  //        ProgressiveDistanceSolverParaboloid<RectangularCuboid> solver(
  //            central_cell, central_moment, 1.0e-14, new_paraboloid);
  //        paraboloid_to_fit =
  //            Paraboloid(datum + solver.getDistance() * new_normal, frame,
  //                       aligned_paraboloid.a(), aligned_paraboloid.b());
  //      }

  // We assume 5 param: 3 translations and 2 curvatures

  double penalty_param = 1.0;
  double lag_param = 0.0;
  for (int n = 0; n < 100; ++n) {
    // Initialize Jacobian
    Eigen::Matrix<double, 8, 27> jacobian_transpose;
    Eigen::Matrix<double, 27, 1> f;
    Eigen::Matrix<double, 27, 27> W;
    W.setZero();
    for (int i = 0; i < 27; i++) W(i, i) = stencil_weights[i];
    Eigen::Matrix<double, 8, 1> x = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        paraboloid_to_fit.getAlignedParaboloid().a(),
        paraboloid_to_fit.getAlignedParaboloid().b()};
    std::cout << "x0 = " << x << std::endl;
    ComputeJacobianAndF(paraboloid_to_fit, jacobian_transpose, f,
                        neighborhood_VF);
    //   std::cout << "Initial Jacobian = " << jacobian_transpose <<
    //   std::endl;
    double constraint = f(13, 0);
    double old_constraint = constraint;
    f(13, 0) *= std::sqrt(penalty_param);
    f(13, 0) += lag_param / (2.0 * std::sqrt(penalty_param));
    for (int i = 0; i < 8; ++i) {
      jacobian_transpose(i, 13) *= std::sqrt(penalty_param);
    }

    // Levenberg-Marquardt algo
    double epsilon_1 = 1.0e-6;
    double epsilon_2 = 1.0e-6;
    double epsilon_3 = 1.0e-6;
    double tau = 1.0e-3;
    int k = 0, k_max = 100;
    double nu = 2.0;
    Eigen::Matrix<double, 8, 8> A =
        jacobian_transpose * W * jacobian_transpose.transpose();
    Eigen::Matrix<double, 8, 1> g = jacobian_transpose * W * f;
    bool found = g.lpNorm<Eigen::Infinity>() < epsilon_1;
    double mu = -1.0e15;
    for (int i = 0; i < 8; ++i) mu = std::max(mu, A(i, i));
    mu *= tau;
    while (!found && k < k_max) {
      k++;
      for (int i = 0; i < 8; ++i) {
        A(i, i) += mu;
      }
      auto h = A.llt().solve(g);
      if (h.norm() < (epsilon_2 * x.norm() + epsilon_2)) {
        found = true;
      } else {
        double F = f.transpose() * W * f;
        Eigen::Matrix<double, 8, 1> x_new = x + h;
        auto& datum = paraboloid_to_fit.getDatum();
        auto& frame = paraboloid_to_fit.getReferenceFrame();
        auto& aligned_paraboloid = paraboloid_to_fit.getAlignedParaboloid();
        Pt new_datum = datum + x_new(0, 0) * frame[0] + x_new(1, 0) * frame[1] +
                       x_new(2, 0) * frame[2];
        UnitQuaternion x_rotation(x_new(3, 0), frame[0]);
        UnitQuaternion y_rotation(x_new(4, 0), frame[1]);
        UnitQuaternion z_rotation(x_new(5, 0), frame[2]);
        auto total_rotation = x_rotation * y_rotation * z_rotation;
        total_rotation.normalize();
        const ReferenceFrame new_frame = total_rotation * frame;
        auto new_paraboloid =
            Paraboloid(new_datum, new_frame, x_new(6, 0), x_new(7, 0));
        // ProgressiveDistanceSolverParaboloid<RectangularCuboid> solver(
        //     neighborhood_VF.getCenterCell(),
        //     neighborhood_VF.getCenterCellStoredMoments(), 1.0e-14,
        //     new_paraboloid);
        // Pt new_corrected_datum =
        //     new_datum + solver.getDistance() * new_frame[2];
        Pt new_corrected_datum = new_datum;
        auto new_corrected_paraboloid = Paraboloid(
            new_corrected_datum, new_frame, x_new(6, 0), x_new(7, 0));
        ComputeJacobianAndF(new_corrected_paraboloid, jacobian_transpose, f,
                            neighborhood_VF);
        constraint = f(13, 0);
        f(13, 0) *= std::sqrt(penalty_param);
        f(13, 0) += lag_param / (2.0 * std::sqrt(penalty_param));
        for (int i = 0; i < 8; ++i) {
          jacobian_transpose(i, 13) *= std::sqrt(penalty_param);
        }
        F -= f.transpose() * W * f;
        double rho = F / (h.transpose() * (mu * h + g));
        if (rho > 0.0) {
          paraboloid_to_fit = Paraboloid(new_corrected_datum, new_frame,
                                         x_new(6, 0), x_new(7, 0));
          x(6, 0) = x_new(6, 0);
          x(7, 0) = x_new(7, 0);
          A = jacobian_transpose * W * jacobian_transpose.transpose();
          g = jacobian_transpose * W * f;
          found = (g.lpNorm<Eigen::Infinity>() < epsilon_1) ||
                  (f.squaredNorm() < epsilon_3);
          mu *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * rho - 1.0, 3.0));
          nu = 2.0;

          //////////////////////////////
          surfaces.clear();
          // const auto cell = IRL::RectangularCuboid::fromBoundingPts(
          //     Pt(-0.5, -0.5, -0.5), Pt(0.5, 0.5, 0.5));
          const auto cell = neighborhood_VF.getCenterCell();
          auto volume_and_surface = getVolumeMoments<
              AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
              cell, new_corrected_paraboloid);
          auto& surface = volume_and_surface.getSurface();
          surface.setLengthScale(std::pow(cell.calculateVolume(), 1.0 / 3.0) /
                                 30.0);
          surfaces.push_back(surface);
          vtk_io.writeVTKInterface(time, surfaces, true);
          ++viz_output;
          /////////////////////////////////
        } else {
          for (int i = 0; i < 8; ++i) {
            A(i, i) -= mu;
          }
          mu *= nu;
          nu *= 2.0;
        }
      }
      //   std::cout << "x = " << x << std::endl;
    }
    std::cout << "Converged in " << k << " iterations" << std::endl;
    std::cout << "Lagrange param = " << lag_param << std::endl;
    std::cout << "Penalty  param = " << penalty_param << std::endl;
    std::cout << "Constraint = " << std::fabs(constraint) << std::endl;

    // Update lagrange multiplier
    lag_param += 2.0 * penalty_param * constraint;
    if (std::fabs(constraint) > 0.25 * std::fabs(old_constraint)) {
      penalty_param *= 2.0;
    }
    old_constraint = constraint;
    if (std::fabs(constraint) < epsilon_1) {
      std::cout << "Exiting big loop after " << n + 1 << " iterations"
                << std::endl;
      break;
    }
  }

  //   // Final result
  //   auto& datum = paraboloid_to_fit.getDatum();
  //   auto& frame = paraboloid_to_fit.getReferenceFrame();
  //   auto& aligned_paraboloid = paraboloid_to_fit.getAlignedParaboloid();
  //   Pt new_datum =
  //       datum + x(0, 0) * frame[0] + x(1, 0) * frame[1] + x(2, 0) *
  //       frame[2];
  //   UnitQuaternion x_rotation(x(3, 0), frame[0]);
  //   UnitQuaternion y_rotation(x(4, 0), frame[1]);
  //   UnitQuaternion z_rotation(x(5, 0), frame[2]);
  //   auto total_rotation = x_rotation * y_rotation * z_rotation;
  //   total_rotation.normalize();
  //   ReferenceFrame new_frame = total_rotation * frame;
  //   auto new_paraboloid = Paraboloid(new_datum, new_frame, x(6, 0), x(7,
  //   0)); surfaces.clear();
  //   // const auto cell = IRL::RectangularCuboid::fromBoundingPts(
  //   //     Pt(-0.5, -0.5, -0.5), Pt(0.5, 0.5, 0.5));
  //   const auto cell = neighborhood_VF.getCenterCell();
  //   auto volume_and_surface =
  //       getVolumeMoments<AddSurfaceOutput<Volume,
  //       ParametrizedSurfaceOutput>>(
  //           cell, new_paraboloid);
  //   auto& surface = volume_and_surface.getSurface();
  //   surface.setLengthScale(std::pow(cell.calculateVolume(), 1.0 / 3.0)
  //   / 30.0); surfaces.push_back(surface); vtk_io.writeVTKInterface(time,
  //   surfaces, true);
  //   ++viz_output;
  //   //    }
}

void ComputeJacobianAndF(Paraboloid& paraboloid,
                         Eigen::Matrix<double, 8, 27>& jacobian_transpose,
                         Eigen::Matrix<double, 27, 1>& f,
                         LVIRANeighborhood<RectangularCuboid> lvira_stencil) {
  using MyGradientType = ParaboloidGradientLocal;
  using MyPtType = PtWithGradient<MyGradientType>;

  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        const int cell_id = i + j * 3 + k * 9;
        auto cell = lvira_stencil.getCell(cell_id);
        auto moments =
            cell.calculateVolume() * lvira_stencil.getStoredMoments(cell_id);
        auto data_cell =
            StoredRectangularCuboid<MyPtType>::fromOtherPolytope(cell);
        auto moments_with_gradients =
            getVolumeMoments<VolumeWithGradient<MyGradientType>>(data_cell,
                                                                 paraboloid);
        f(cell_id, 0) = moments - moments_with_gradients.volume();
        auto gradient = moments_with_gradients.volume_gradient();
        jacobian_transpose(0, cell_id) = gradient.getGradTx();
        jacobian_transpose(1, cell_id) = gradient.getGradTy();
        jacobian_transpose(2, cell_id) = gradient.getGradTz();
        jacobian_transpose(3, cell_id) = gradient.getGradRx();
        jacobian_transpose(4, cell_id) = gradient.getGradRy();
        jacobian_transpose(5, cell_id) = gradient.getGradRz();
        jacobian_transpose(6, cell_id) = gradient.getGradA();
        jacobian_transpose(7, cell_id) = gradient.getGradB();
      }
    }
  }
}

}  // namespace IRL
