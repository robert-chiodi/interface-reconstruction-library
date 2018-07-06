// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_ADVECTED_PLANE_RECONSTRUCTION_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_ADVECTED_PLANE_RECONSTRUCTION_TPP_

namespace IRL {

template <class MomentsContainerType, class CellType>
PlanarSeparator AdvectedPlaneReconstruction::solve(
    const MomentsContainerType& a_volume_moments_list,
    const R2PNeighborhood<CellType>& a_neighborhood, const double threshold) {
  Normal mean_normal = Normal(0.0, 0.0, 0.0);
  double summed_volume = {0.0};
  for (const auto& element : a_volume_moments_list) {
    mean_normal += element.normal();
    summed_volume += element.volumeMoments().volume();
  }
  mean_normal /= summed_volume;
  double mean_vector_magnitude = magnitude(mean_normal);
  double target_volume_fraction =
      a_neighborhood.getCenterCellStoredMoments()[0].volume() /
      a_neighborhood.getCenterCell().calculateVolume();

  // If magnitude of this vector > threshold, one plane is sufficient
  if (mean_vector_magnitude > threshold) {
    mean_normal /= mean_vector_magnitude;
    return PlanarSeparator::fromOnePlane(
        Plane(mean_normal,
              findDistanceOnePlane(a_neighborhood.getCenterCell(),
                                   target_volume_fraction, mean_normal)));
  }
  // If we made it this far, need atleast two planes
  // to represent listed moments.

  // Partition moments according to normal vectors.
  PartitionByNormal<MomentsContainerType> moment_partitioning(
      &a_volume_moments_list);
  moment_partitioning.setup();
  KMeans::partition(&moment_partitioning);
  auto partitioned_moments = moment_partitioning.getPartitionedObjects();

  // Test four permutations of found normal/centroids. Accept one with lowest
  // error on the R2PNeighborhood.
  return AdvectedPlaneReconstruction::findBestPermutation(
      &partitioned_moments, a_neighborhood, target_volume_fraction);
}

template <class ContainedType, class CellType>
PlanarSeparator AdvectedPlaneReconstruction::findBestPermutation(
    ContainedType* a_moments, const R2PNeighborhood<CellType>& a_neighborhood,
    const double a_target_volume_fraction) {
  // Normalize centroids and normals.
  (*a_moments)[0].normalizeByVolume();
  (*a_moments)[0].normal().normalize();
  (*a_moments)[1].normalizeByVolume();
  (*a_moments)[1].normal().normalize();

  PlanarSeparator best_separator, attempted_separator;
  double min_err = DBL_MAX;
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[0].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[1].volumeMoments().centroid(),
      1.0);
  AdvectedPlaneReconstruction::checkIfBest(a_neighborhood, attempted_separator,
                                           &best_separator, &min_err);
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[0].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[1].volumeMoments().centroid(),
      -1.0);
  AdvectedPlaneReconstruction::checkIfBest(a_neighborhood, attempted_separator,
                                           &best_separator, &min_err);
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[1].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[0].volumeMoments().centroid(),
      1.0);
  AdvectedPlaneReconstruction::checkIfBest(a_neighborhood, attempted_separator,
                                           &best_separator, &min_err);
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[1].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[0].volumeMoments().centroid(),
      -1.0);
  AdvectedPlaneReconstruction::checkIfBest(a_neighborhood, attempted_separator,
                                           &best_separator, &min_err);
  return best_separator;
}

template <class CellType>
PlanarSeparator AdvectedPlaneReconstruction::constructSeparatorAttempt(
    const CellType& a_cell, const double a_target_volume_fraction,
    const Normal& a_normal_0, const Pt& a_pt_0, const Normal& a_normal_1,
    const Pt& a_pt_1, const double a_flip_cut) {
  // Construct PlanarSeparator
  PlanarSeparator attempted_separator = PlanarSeparator::fromTwoPlanes(
      Plane(a_normal_0, a_normal_0 * a_pt_0),
      Plane(a_normal_1, a_normal_1 * a_pt_1), a_flip_cut);
  // Find volume conserving distance
  IterativeSolverForDistance<CellType> solver(
      a_cell, a_target_volume_fraction,
      global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE,
      attempted_separator);
  attempted_separator.setDistances(solver.getDistances());
  // Clean the reconstruction
  cleanReconstruction(a_cell, a_target_volume_fraction, &attempted_separator);
  return attempted_separator;
}

template <class CellType>
void AdvectedPlaneReconstruction::checkIfBest(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_attempt_separator,
    PlanarSeparator* a_current_best_separator,
    double* a_current_minimum_error) {
  double err = 0.0;
  // Iterate over CellGroupedMoments in a_neighborhood and accumulate
  // error in centroids.
  UnsignedIndex_t ind = 0;
  for (const auto& cell_grouped_moments : a_neighborhood) {
    auto svm = cell_grouped_moments.calculateNormalizedVolumeMoments(
        a_attempt_separator);
    auto cell_volume = cell_grouped_moments.getStoredMoments()[0].volume() +
                       cell_grouped_moments.getStoredMoments()[1].volume();
    auto cell_VF = svm[0].volume() / cell_volume;
    if (cell_VF < global_constants::VF_LOW) {
      svm[0].centroid() =
          a_neighborhood.getCenterCellStoredMoments()[0].centroid();
    }
    if (cell_VF > global_constants::VF_HIGH) {
      svm[1].centroid() =
          a_neighborhood.getCenterCellStoredMoments()[1].centroid();
    }

    if (cell_grouped_moments.getStoredMoments()[0].volume() / cell_volume >
        global_constants::VF_LOW) {
      err += magnitude(cell_grouped_moments.getStoredMoments()[0].centroid() -
                       svm[0].centroid());  // Liquid centroid contribution
    }
    if (cell_grouped_moments.getStoredMoments()[1].volume() / cell_volume >
        global_constants::VF_LOW) {
      err += magnitude(cell_grouped_moments.getStoredMoments()[1].centroid() -
                       svm[1].centroid());  // Gas centroid contribution
    }
    ++ind;
  }
  if (err < *a_current_minimum_error) {
    *a_current_minimum_error = err;
    *a_current_best_separator = a_attempt_separator;
  }
}

//******************************************************************* //
//     Debugging function template implementations below this
//******************************************************************* //

template <class MomentsContainerType, class CellType>
PlanarSeparator AdvectedPlaneReconstructionDebug::solve(
    const MomentsContainerType& a_volume_moments_list,
    const R2PNeighborhood<CellType>& a_neighborhood, const double threshold) {
  AdvectedPlaneReconstructionDebug::writeOutCentroids(a_neighborhood);
  Normal mean_normal = Normal(0.0, 0.0, 0.0);
  double summed_volume = {0.0};
  for (const auto& element : a_volume_moments_list) {
    mean_normal += element.normal();
    summed_volume += element.volumeMoments().volume();
  }
  mean_normal /= summed_volume;
  double mean_vector_magnitude = magnitude(mean_normal);
  double target_volume_fraction =
      a_neighborhood.getCenterCellStoredMoments()[0].volume() /
      a_neighborhood.getCenterCell().calculateVolume();
  // If magnitude of this vector > threshold, one plane is sufficient
  if (mean_vector_magnitude > threshold) {
    std::cout << "% Mean vector magnitude : " << mean_vector_magnitude << '\n';
    std::cout << "% Total surface density : "
              << summed_volume /
                     a_neighborhood.getCenterCell().calculateVolume()
              << '\n';
    mean_normal.normalize();
    PlanarSeparator reconstruction_to_return = PlanarSeparator::fromOnePlane(
        Plane(mean_normal,
              findDistanceOnePlane(a_neighborhood.getCenterCell(),
                                   target_volume_fraction, mean_normal)));
    AdvectedPlaneReconstructionDebug::writeOutPlane(
        a_neighborhood, reconstruction_to_return, "bestPlane", 0);
    return reconstruction_to_return;
  }
  // If we made it this far, need atleast two planes
  // to represent listed moments.

  // Partition moments according to normal vectors.
  PartitionByNormal<MomentsContainerType> moment_partitioning(
      &a_volume_moments_list);
  moment_partitioning.setup();
  KMeans::partition(&moment_partitioning);
  auto partitioned_moments = moment_partitioning.getPartitionedObjects();

  // Test four permutations of found normal/centroids. Accept one with lowest
  // error on the R2PNeighborhood.
  return AdvectedPlaneReconstructionDebug::findBestPermutation(
      &partitioned_moments, a_neighborhood, target_volume_fraction);
}

template <class ContainedType, class CellType>
PlanarSeparator AdvectedPlaneReconstructionDebug::findBestPermutation(
    ContainedType* a_moments, const R2PNeighborhood<CellType>& a_neighborhood,
    const double a_target_volume_fraction) {
  // Normalize centroids and normals.
  (*a_moments)[0].normalizeByVolume();
  (*a_moments)[0].normal().normalize();
  (*a_moments)[1].normalizeByVolume();
  (*a_moments)[1].normal().normalize();

  PlanarSeparator best_separator, attempted_separator;
  double min_err = DBL_MAX;
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[0].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[1].volumeMoments().centroid(),
      1.0);
  AdvectedPlaneReconstructionDebug::writeOutPlane(
      a_neighborhood, attempted_separator, "bestPlane", 0);
  AdvectedPlaneReconstructionDebug::checkIfBest(
      a_neighborhood, attempted_separator, &best_separator, &min_err);
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[0].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[1].volumeMoments().centroid(),
      -1.0);
  AdvectedPlaneReconstructionDebug::writeOutPlane(
      a_neighborhood, attempted_separator, "bestPlane", 1);
  AdvectedPlaneReconstructionDebug::checkIfBest(
      a_neighborhood, attempted_separator, &best_separator, &min_err);
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[1].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[0].volumeMoments().centroid(),
      1.0);
  AdvectedPlaneReconstructionDebug::writeOutPlane(
      a_neighborhood, attempted_separator, "bestPlane", 2);
  AdvectedPlaneReconstructionDebug::checkIfBest(
      a_neighborhood, attempted_separator, &best_separator, &min_err);
  attempted_separator = AdvectedPlaneReconstruction::constructSeparatorAttempt(
      a_neighborhood.getCenterCell(), a_target_volume_fraction,
      (*a_moments)[0].normal(), (*a_moments)[1].volumeMoments().centroid(),
      (*a_moments)[1].normal(), (*a_moments)[0].volumeMoments().centroid(),
      -1.0);
  AdvectedPlaneReconstructionDebug::writeOutPlane(
      a_neighborhood, attempted_separator, "bestPlane", 3);
  AdvectedPlaneReconstructionDebug::checkIfBest(
      a_neighborhood, attempted_separator, &best_separator, &min_err);
  return best_separator;
}

template <class CellType>
void AdvectedPlaneReconstructionDebug::checkIfBest(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_attempt_separator,
    PlanarSeparator* a_current_best_separator,
    double* a_current_minimum_error) {
  double err = 0.0;
  // Iterate over CellGroupedMoments in a_neighborhood and accumulate
  // error in centroids.
  UnsignedIndex_t ind = 0;
  for (const auto& cell_grouped_moments : a_neighborhood) {
    auto svm = cell_grouped_moments.calculateNormalizedVolumeMoments(
        a_attempt_separator);
    auto cell_volume = cell_grouped_moments.getStoredMoments()[0].volume() +
                       cell_grouped_moments.getStoredMoments()[1].volume();
    auto cell_VF = svm[0].volume() / cell_volume;
    if (cell_VF < global_constants::VF_LOW) {
      svm[0].centroid() =
          a_neighborhood.getCenterCellStoredMoments()[0].centroid();
    }
    if (cell_VF > global_constants::VF_HIGH) {
      svm[1].centroid() =
          a_neighborhood.getCenterCellStoredMoments()[1].centroid();
    }

    if (cell_grouped_moments.getStoredMoments()[0].volume() / cell_volume >
        global_constants::VF_LOW) {
      err += magnitude(cell_grouped_moments.getStoredMoments()[0].centroid() -
                       svm[0].centroid());  // Liquid centroid contribution
    }
    if (cell_grouped_moments.getStoredMoments()[1].volume() / cell_volume >
        global_constants::VF_LOW) {
      err += magnitude(cell_grouped_moments.getStoredMoments()[1].centroid() -
                       svm[1].centroid());  // Gas centroid contribution
    }
    ++ind;
  }
  std::cout << "% Error for reconstruction: " << err << '\n';
  if (err < *a_current_minimum_error) {
    *a_current_minimum_error = err;
    *a_current_best_separator = a_attempt_separator;
  }
}  // namespace IRL

template <class CellType>
void AdvectedPlaneReconstructionDebug::writeOutPlane(
    const R2PNeighborhood<CellType>& a_neighborhood,
    const PlanarSeparator& a_reconstruction, const std::string& a_prefix,
    const std::size_t a_iteration_number) {
  PlanarSeparator tmp = a_reconstruction;
  std::cout << "%**************************************** \n";
  std::cout << "% Writing out reconstruction: " << a_iteration_number << '\n';
  std::cout << "% Saved reconstruction has " << tmp.getNumberOfPlanes()
            << "planes \n";
  for (UnsignedIndex_t n = 0; n < tmp.getNumberOfPlanes(); ++n) {
    std::cout << "% Plane " << n << ":  " << tmp[n].normal()[0] << " "
              << tmp[n].normal()[1] << " " << tmp[n].normal()[2] << " "
              << tmp[n].distance() << '\n';
  }
  std::cout << "% Is reconstruction flipped?  " << tmp.isFlipped() << '\n';
  std::cout << "% **************************************** \n";
  for (UnsignedIndex_t cell = 0; cell < a_neighborhood.size(); ++cell) {
    for (UnsignedIndex_t plane = 0;
         plane < a_reconstruction.getNumberOfPlanes(); ++plane) {
      Polygon poly = getPlanePolygonFromReconstruction<Polygon>(
          a_neighborhood.getCell(cell), a_reconstruction,
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

template <class CellType>
void AdvectedPlaneReconstructionDebug::writeOutCentroids(
    const R2PNeighborhood<CellType>& a_neighborhood) {
  Pt center_cell_centroid = a_neighborhood.getCenterCell().calculateCentroid();
  std::cout << "x0 = " << center_cell_centroid[0] << ";\n";
  std::cout << "y0 = " << center_cell_centroid[1] << ";\n";
  std::cout << "z0 = " << center_cell_centroid[2] << ";\n";
  for (UnsignedIndex_t n = 0; n < a_neighborhood.size(); ++n) {
    std::cout << "liq_centroid(1:4," << n + 1 << ") = [";
    Pt liq_centroid(a_neighborhood.getStoredMoments(n)[0].centroid()[0],
                    a_neighborhood.getStoredMoments(n)[0].centroid()[1],
                    a_neighborhood.getStoredMoments(n)[0].centroid()[2]);
    std::cout << liq_centroid.x() << "," << liq_centroid.y() << ","
              << liq_centroid.z() << "," << 1.0 << "];\n";
    std::cout << "gas_centroid(1:4," << n + 1 << ") = [";
    Pt gas_centroid(a_neighborhood.getStoredMoments(n)[1].centroid()[0],
                    a_neighborhood.getStoredMoments(n)[1].centroid()[1],
                    a_neighborhood.getStoredMoments(n)[1].centroid()[2]);
    std::cout << gas_centroid.x() << "," << gas_centroid.y() << ","
              << gas_centroid.z() << "," << 1.0 << "];\n";
  }
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_ADVECTED_PLANE_RECONSTRUCTION_TPP_
