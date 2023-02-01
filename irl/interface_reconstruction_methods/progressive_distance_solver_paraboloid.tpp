// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_TPP_

#include <algorithm>

#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"

namespace IRL {

template <class CellType>
ProgressiveDistanceSolverParaboloid<CellType>::
    ProgressiveDistanceSolverParaboloid(
        const CellType& a_cell, const double a_volume_fraction,
        const double a_volume_fraction_tolerance,
        const Paraboloid& a_reconstruction)
    : target_volume_fraction_m(a_volume_fraction),
      volume_fraction_tolerance_m(a_volume_fraction_tolerance),
      reconstruction_m(a_reconstruction) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  this->solveForDistance(a_cell);
}

template <class CellType>
void ProgressiveDistanceSolverParaboloid<CellType>::solve(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const Paraboloid& a_reconstruction) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  target_volume_fraction_m = a_volume_fraction;
  volume_fraction_tolerance_m = a_volume_fraction_tolerance;
  reconstruction_m = a_reconstruction;
  this->solveForDistance(a_cell);
}

template <class CellType>
double ProgressiveDistanceSolverParaboloid<CellType>::getDistance(void) {
  return distances_m;
}

template <class CellType>
bool ProgressiveDistanceSolverParaboloid<CellType>::isBoundsTrueBounds(
    const CellType& a_cell) {
  const auto datum = reconstruction_m.getDatum();
  const auto& ref_frame = reconstruction_m.getReferenceFrame();
  const auto& aligned_paraboloid = reconstruction_m.getAlignedParaboloid();

  const Pt low_datum = datum + sorted_distances_m.front() * ref_frame[2];
  std::cout << "Low datum = " << low_datum << std::endl;
  const Pt high_datum = datum + sorted_distances_m.back() * ref_frame[2];
  std::cout << "High datum = " << high_datum << std::endl;
  const auto low_paraboloid = Paraboloid(
      low_datum, ref_frame, aligned_paraboloid.a(), aligned_paraboloid.b());
  const auto high_paraboloid = Paraboloid(
      high_datum, ref_frame, aligned_paraboloid.a(), aligned_paraboloid.b());
  double result = getVolumeFraction(a_cell, low_paraboloid);
  if (result > target_volume_fraction_m) {
    // Volume fraction should be less than target_volume_fraction_m
    std::cout << "Lower bound result failed: " << result << std::endl;
    return false;
  } else {
    std::cout << "Lower bound result passed: " << result << std::endl;
  }
  result = getVolumeFraction(a_cell, high_paraboloid);
  if (result < target_volume_fraction_m) {
    // Volume fraction should be > target_volume_fraction_m.
    std::cout << "Upper bound result failed: " << result << std::endl;
    return false;
  } else {
    std::cout << "Upper bound result passed: " << result << std::endl;
  }
  return true;
}

// template <class CellType>
// void ProgressiveDistanceSolverParaboloid<CellType>::solveForDistance(
//     const CellType& a_cell) {
//   // Calculate volume of cell
//   auto initial_cell_volume = a_cell.calculateVolume();

//   // Generate half-edge structure
//   auto& complete_polytope = setHalfEdgeStructureParaboloid(a_cell);
//   auto half_edge_polytope =
//       generateSegmentedVersionParaboloid<CellType>(&complete_polytope);
//   assert(half_edge_polytope.checkValidHalfEdgeStructure());
//   const UnsignedIndex_t number_of_vertices =
//       half_edge_polytope.getNumberOfVertices();
//   const double inv_scale =
//       std::max(1.0e-6, std::pow(initial_cell_volume, 1.0 / 3.0));
//   const double scale = 1.0 / inv_scale;

//   // Move cell to local frame of reference of the paraboloid
//   const auto& datum = scale * reconstruction_m.getDatum();
//   const auto& ref_frame = reconstruction_m.getReferenceFrame();
//   const auto aligned_paraboloid = AlignedParaboloid(std::array<double, 2>(
//       {inv_scale * reconstruction_m.getAlignedParaboloid().a(),
//        inv_scale * reconstruction_m.getAlignedParaboloid().b()}));
//   for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
//     const Pt original_pt =
//         scale * half_edge_polytope.getVertex(v)->getLocation().getPt() -
//         datum;
//     typename CellType::pt_type projected_location;
//     auto& pt = projected_location.getPt();
//     for (UnsignedIndex_t n = 0; n < 3; ++n) {
//       pt[n] = ref_frame[n] * original_pt;
//     }
//     // std::cout << "New pt = " << pt << std::endl;
//     half_edge_polytope.getVertex(v)->setLocation(projected_location);
//   }

//   initial_cell_volume_m = half_edge_polytope.calculateVolume();

//   // Sort the vertices Z positions in the local frame of reference
//   std::vector<Pt> original_vertices;
//   original_vertices.resize(number_of_vertices);
//   for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
//     original_vertices[v] =
//         half_edge_polytope.getVertex(v)->getLocation().getPt();
//   }

//   // Compute new face planes
//   for (auto& face : half_edge_polytope) {
//     auto normal = Normal(0.0, 0.0, 0.0);
//     const auto starting_half_edge = face->getStartingHalfEdge();
//     auto current_half_edge = starting_half_edge;
//     auto next_half_edge = starting_half_edge->getNextHalfEdge();
//     const auto& start_location =
//         starting_half_edge->getPreviousVertex()->getLocation().getPt();
//     do {
//       normal += crossProduct(
//           current_half_edge->getVertex()->getLocation().getPt() -
//               start_location,
//           next_half_edge->getVertex()->getLocation().getPt() -
//           start_location);
//       current_half_edge = next_half_edge;
//       next_half_edge = next_half_edge->getNextHalfEdge();
//     } while (next_half_edge != starting_half_edge);
//     normal.normalize();
//     face->setPlane(Plane(normal, normal * start_location));
//   }
//   assert(half_edge_polytope.checkValidHalfEdgeStructure());

//   // Compute and sort distances between the paraboloid and the vertices
//   // ---> in the local frame of reference of the paraboloid
//   sorted_distances_m.resize(number_of_vertices + 2);
//   UnsignedIndex_t count = 0;
//   const Pt& start_pt =
//   half_edge_polytope.getVertex(0)->getLocation().getPt();
//   sorted_distances_m[count++] = signedDistance(start_pt, aligned_paraboloid);
//   double cell_zmin = start_pt[2], cell_zmax = start_pt[2];
//   for (UnsignedIndex_t v = 1; v < number_of_vertices; ++v) {
//     const Pt& vertex =
//     half_edge_polytope.getVertex(v)->getLocation().getPt();
//     sorted_distances_m[count++] = signedDistance(vertex, aligned_paraboloid);
//     cell_zmin = std::min({vertex[2], cell_zmin});
//     cell_zmax = std::max({vertex[2], cell_zmax});
//   }
//   // Sort  into ascending
//   sorted_distances_m[count++] = -1.0e15;
//   sorted_distances_m[count++] = 1.0e15;
//   std::sort(sorted_distances_m.begin(), sorted_distances_m.end());
//   sorted_distances_m[0] =
//       sorted_distances_m[1] - 10.0 * (cell_zmax - cell_zmin);
//   sorted_distances_m[number_of_vertices + 1] =
//       sorted_distances_m[number_of_vertices] + 10.0 * (cell_zmax -
//       cell_zmin);

//   // std::cout << "List of sorted distances = ";
//   // for (UnsignedIndex_t v = 0; v < number_of_vertices + 2; ++v) {
//   //   std::cout << sorted_distances_m[v] << "  ";
//   // }
//   // std::cout << std::endl;
//   // assert(this->isBoundsTrueBounds(a_cell));

//   // Keep bisectioning between nodes until we find between where the solution
//   // lays. Will then just have a prismatoid that needs to be optimized over
//   std::array<UnsignedIndex_t, 3> bounding_indices{
//       {0, 0, static_cast<UnsignedIndex_t>(sorted_distances_m.size()) - 1}};
//   std::array<double, 3> bounding_values{{0.0, 0.0, 1.0}};
//   // while (bounding_indices[2] - bounding_indices[0] > 1) {
//   //   bounding_indices[1] = (bounding_indices[0] + bounding_indices[2]) / 2;
//   //   // std::cout << "Distance = " <<
//   sorted_distances_m[bounding_indices[1]]
//   //   // << std::endl;
//   //   shiftSegmentedPolytopeAndUpdateFaces(
//   //       &half_edge_polytope, original_vertices,
//   //       sorted_distances_m[bounding_indices[1]]);
//   //   assert(half_edge_polytope.checkValidHalfEdgeStructure());
//   //   // std::cout << "Caculating volume... " << std::endl;
//   //   const double volume_fraction_cut =
//   //       intersectPolyhedronWithParaboloid<Volume>(
//   //           &half_edge_polytope, &complete_polytope, aligned_paraboloid) /
//   //       initial_cell_volume_m;
//   //   if (volume_fraction_cut * initial_cell_volume_m < -999.0) {
//   //     distances_m =
//   //         -999999.0 + inv_scale * sorted_distances_m[bounding_indices[1]];
//   //     std::cout << "Initial bisection between " << sorted_distances_m[0]
//   //               << " and " << sorted_distances_m[number_of_vertices + 1]
//   //               << " / distance = "
//   //               << inv_scale * sorted_distances_m[bounding_indices[1]]
//   //               << std::endl;
//   //     return;
//   //   }
//   //   resetPolyhedron(&half_edge_polytope, &complete_polytope);
//   //   // std::cout << "Volume calculated is " << volume_fraction_cut <<
//   //   std::endl; if (volume_fraction_cut >
//   //       target_volume_fraction_m + volume_fraction_tolerance_m) {
//   //     bounding_indices[2] = bounding_indices[1];
//   //     bounding_values[2] = volume_fraction_cut;
//   //   } else if (volume_fraction_cut <
//   //              target_volume_fraction_m - volume_fraction_tolerance_m) {
//   //     bounding_indices[0] = bounding_indices[1];
//   //     bounding_values[0] = volume_fraction_cut;
//   //   } else {
//   //     // std::cout << "DISTANCE1 = " <<
//   //     sorted_distances_m[bounding_indices[1]]; distances_m = inv_scale *
//   //     sorted_distances_m[bounding_indices[1]]; return;
//   //   }
//   // }

//   double interval_min = sorted_distances_m[bounding_indices[0]];
//   double interval_max = sorted_distances_m[bounding_indices[2]];

//   // std::cout << std::endl
//   //           << "Searching in the interval : ["
//   //           << sorted_distances_m[bounding_indices[0]] << ", "
//   //           << sorted_distances_m[bounding_indices[2]] << "]" << std::endl
//   //           << std::endl;

//   // // Now just have a prismatoid left. Perform secant, fall back to
//   // Bisection if it fails.
//   // double distance = sorted_distances_m[bounding_indices[0]];
//   // double delta = sorted_distances_m[bounding_indices[0]] -
//   //                sorted_distances_m[bounding_indices[2]];
//   // double old_error = bounding_values[2] - target_volume_fraction_m;
//   // double error = bounding_values[0] - target_volume_fraction_m;
//   // bounding_values[0] = sorted_distances_m[bounding_indices[0]];
//   // bounding_values[2] = sorted_distances_m[bounding_indices[2]];
//   // for (UnsignedIndex_t iter = 0; iter < max_iter_m; ++iter) {
//   //   delta *= -error / safelyEpsilon(error - old_error);
//   //   old_error = error;
//   //   distance += delta;
//   //   if (distance < interval_min || distance > interval_max) {
//   //     break;
//   //   }
//   //   // std::cout << "Distance = " << distance << std::endl;
//   //   shiftSegmentedPolytopeAndUpdateFaces(&half_edge_polytope,
//   //   original_vertices,
//   //                                        distance);
//   //   // std::cout << "Caculating volume... " << std::endl;
//   //   const double volume_fraction_cut =
//   //       intersectPolyhedronWithParaboloid<Volume>(
//   //           &half_edge_polytope, &complete_polytope, aligned_paraboloid) /
//   //       initial_cell_volume_m;
//   //   if (volume_fraction_cut * initial_cell_volume_m < -999.0) {
//   //     distances_m = -999999.0 + inv_scale * distance;
//   //     std::cout << "Secant between " << interval_min << " and " <<
//   //     interval_max
//   //               << " / distance = " << inv_scale * distance << std::endl;
//   //     return;
//   //   }
//   //   resetPolyhedron(&half_edge_polytope, &complete_polytope);
//   //   // std::cout << "Volume calculated is " << volume_fraction_cut
//   //   //           << " for target " << target_volume_fraction_m <<
//   std::endl;
//   //   error = volume_fraction_cut - target_volume_fraction_m;
//   //   if (volume_fraction_cut >
//   //       target_volume_fraction_m + volume_fraction_tolerance_m) {
//   //     if (distance < bounding_values[2]) {
//   //       bounding_values[2] = distance;
//   //     }
//   //   } else if (volume_fraction_cut <
//   //              target_volume_fraction_m - volume_fraction_tolerance_m) {
//   //     if (distance > bounding_values[0]) {
//   //       bounding_values[0] = distance;
//   //     }
//   //   } else if (!std::isnan(volume_fraction_cut)) {
//   //     distances_m = inv_scale * distance;
//   //     return;
//   //   }
//   // }
//   bounding_values[0] = interval_min;  // std::max(bounding_values[0],
//                                       // sorted_distances_m.front());
//   bounding_values[2] =
//       interval_max;  // std::min(bounding_values[2],
//       sorted_distances_m.back());

//   // std::cout << std::endl
//   //           << "Performing bisection (fall-back)in the interval : ["
//   //           << bounding_values[0] << ", " << bounding_values[2] << "]"
//   //           << std::endl
//   //           << std::endl;

//   // Perform bisection since secant failed to find answer within tolerance.
//   // Move cell back to its initial position

//   for (UnsignedIndex_t iter = 0; iter < max_bisection_iter; ++iter) {
//     bounding_values[1] = 0.5 * (bounding_values[0] + bounding_values[2]);
//     // std::cout << "Distance =  " << bounding_values[1] << std::endl;
//     shiftSegmentedPolytopeAndUpdateFaces(&half_edge_polytope,
//     original_vertices,
//                                          bounding_values[1]);
//     // std::cout << "Caculating volume... " << std::endl;
//     const double volume_fraction_cut =
//         intersectPolyhedronWithAlignedParaboloid<Volume>(
//             &half_edge_polytope, &complete_polytope, aligned_paraboloid) /
//         initial_cell_volume_m;
//     if (volume_fraction_cut * initial_cell_volume_m < -999.0) {
//       distances_m = -999999.0 + inv_scale * bounding_values[1];
//       std::cout << "Bisection between " << bounding_values[0] << " and "
//                 << bounding_values[2]
//                 << " / distance = " << inv_scale * bounding_values[1]
//                 << std::endl;
//       return;
//     }
//     resetPolyhedron(&half_edge_polytope, &complete_polytope);
//     // std::cout << "Volume calculated is " << volume_fraction_cut <<
//     std::endl; if (volume_fraction_cut >
//         target_volume_fraction_m + volume_fraction_tolerance_m) {
//       bounding_values[2] = bounding_values[1];
//     } else if (volume_fraction_cut <
//                target_volume_fraction_m - volume_fraction_tolerance_m) {
//       bounding_values[0] = bounding_values[1];
//     } else {
//       // std::cout << "DISTANCE3 = " << bounding_values[1];
//       distances_m = inv_scale * bounding_values[1];
//       return;
//     }
//   }

//   const double volume_fraction_cut =
//       intersectPolyhedronWithAlignedParaboloid<Volume>(
//           &half_edge_polytope, &complete_polytope, aligned_paraboloid) /
//       initial_cell_volume_m;
//   if (volume_fraction_cut * initial_cell_volume_m < -999.0) {
//     std::string poly_filename = "error_cell";
//     std::cout << "Bisection FAILED between " << bounding_values[0] << " and "
//               << bounding_values[2]
//               << " / distance = " << inv_scale * distances_m << std::endl;
//     distances_m = -999999.0 + inv_scale * distances_m;
//     return;
//   }

//   // if (std::fabs(volume_fraction_cut - target_volume_fraction_m) > 1.0e-13)
//   {
//   //   std::cout << "BISECTION HAS FAILED : volume fraction = "
//   //             << volume_fraction_cut << " instead of "
//   //             << target_volume_fraction_m
//   //             << " diff = " << (volume_fraction_cut -
//   //             target_volume_fraction_m)
//   //             << "\n Zmin/max = " << sorted_distances_m[0] << " -- "
//   //             << sorted_distances_m[number_of_vertices + 1]
//   //             << "\n Bounding values = " << bounding_values[0] << " -- "
//   //             << bounding_values[1] << " -- " << bounding_values[2]
//   //             << "\n Paraboloid = " << aligned_paraboloid.a() << " -- "
//   //             << aligned_paraboloid.b() << std::endl;
//   // }
//   std::cout << "DISTANCE FAILED = " << inv_scale * bounding_values[1];
//   distances_m = -999999.0 + inv_scale * bounding_values[1];
//   return;
// }  // namespace IRL

// template <class SegmentedHalfEdgePolytopeType>
// void shiftSegmentedPolytopeAndUpdateFaces(
//     SegmentedHalfEdgePolytopeType* a_seg_half_edge,
//     const std::vector<Pt>& a_original_verts, const double a_distance) {
//   // Shift along Z
//   const UnsignedIndex_t number_of_vertices =
//       a_seg_half_edge->getNumberOfVertices();
//   // assert(a_original_vert_Z.size() == number_of_vertices);
//   for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
//     // const Pt pt = a_seg_half_edge->getVertex(v)->getLocation().getPt();
//     typename SegmentedHalfEdgePolytopeType::pt_type new_location;
//     Pt& new_pt = new_location.getPt();
//     new_pt[0] = a_original_verts[v][0];
//     new_pt[1] = a_original_verts[v][1];
//     new_pt[2] = a_original_verts[v][2] - a_distance;
//     // std::cout << "New point = " << new_pt << std::endl;
//     a_seg_half_edge->getVertex(v)->setLocation(new_location);
//   }
//   // Update faces
//   for (auto& face : (*a_seg_half_edge)) {
//     const auto normal = face->getPlane().normal();
//     const auto starting_half_edge = face->getStartingHalfEdge();
//     const auto& start_location =
//         starting_half_edge->getPreviousVertex()->getLocation().getPt();
//     // std::cout << "Pt ref = " << start_location << std::endl;
//     // std::cout << "New distance = " << normal * start_location <<
//     std::endl; face->setPlane(Plane(normal, normal * start_location));
//   }
// }

template <class CellType>
void ProgressiveDistanceSolverParaboloid<CellType>::solveForDistance(
    const CellType& a_cell) {
  auto copy_cell = CellType(a_cell);

  // Calculate volume of cell
  auto cell_volume = copy_cell.calculateVolume();
  double length_scale = std::cbrt(cell_volume);

  // Move cell to local frame of reference of the paraboloid
  const auto datum = reconstruction_m.getDatum();
  const auto frame = reconstruction_m.getReferenceFrame();

  double interval_min = -length_scale;
  double interval_max = length_scale;

  reconstruction_m.setDatum(Pt(datum + frame[2] * interval_max));
  double vfrac_max =
      getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
  UnsignedIndex_t iter = 0;
  while (iter < 40 && vfrac_max < (1.0 - volume_fraction_tolerance_m)) {
    interval_max *= 2.0;
    reconstruction_m.setDatum(Pt(datum + frame[2] * interval_max));
    vfrac_max =
        getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
    iter++;
  }

  reconstruction_m.setDatum(Pt(datum + frame[2] * interval_min));
  double vfrac_min =
      getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
  iter = 0;
  while (iter < 40 && vfrac_min > volume_fraction_tolerance_m) {
    interval_min *= 2.0;
    reconstruction_m.setDatum(Pt(datum + frame[2] * interval_min));
    vfrac_min =
        getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
    iter++;
  }

  // Perform bisection since secant failed to find answer within tolerance.
  // Move cell back to its initial position
  std::array<double, 3> bounding_values{{interval_min, 0.0, interval_max}};
  for (UnsignedIndex_t iter = 0; iter < max_bisection_iter; ++iter) {
    bounding_values[1] = 0.5 * (bounding_values[0] + bounding_values[2]);
    reconstruction_m.setDatum(Pt(datum + frame[2] * bounding_values[1]));
    const double vfrac_cut =
        getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
    if (vfrac_cut > target_volume_fraction_m + volume_fraction_tolerance_m) {
      bounding_values[2] = bounding_values[1];
    } else if (vfrac_cut <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      bounding_values[0] = bounding_values[1];
    } else {
      distances_m = bounding_values[1];
      return;
    }
  }

  std::cout << "DISTANCE FAILED = " << bounding_values[1];
  distances_m = -999999.0 + bounding_values[1];
  return;
}  // namespace IRL

template <class SegmentedHalfEdgePolytopeType>
void shiftSegmentedPolytopeAndUpdateFaces(
    SegmentedHalfEdgePolytopeType* a_seg_half_edge,
    const std::vector<Pt>& a_original_verts, const double a_distance) {
  // Shift along Z
  const UnsignedIndex_t number_of_vertices =
      a_seg_half_edge->getNumberOfVertices();
  // assert(a_original_vert_Z.size() == number_of_vertices);
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    // const Pt pt = a_seg_half_edge->getVertex(v)->getLocation().getPt();
    typename SegmentedHalfEdgePolytopeType::pt_type new_location;
    Pt& new_pt = new_location.getPt();
    new_pt[0] = a_original_verts[v][0];
    new_pt[1] = a_original_verts[v][1];
    new_pt[2] = a_original_verts[v][2] - a_distance;
    // std::cout << "New point = " << new_pt << std::endl;
    a_seg_half_edge->getVertex(v)->setLocation(new_location);
  }
  // Update faces
  for (auto& face : (*a_seg_half_edge)) {
    const auto normal = face->getPlane().normal();
    const auto starting_half_edge = face->getStartingHalfEdge();
    const auto& start_location =
        starting_half_edge->getPreviousVertex()->getLocation().getPt();
    // std::cout << "Pt ref = " << start_location << std::endl;
    // std::cout << "New distance = " << normal * start_location << std::endl;
    face->setPlane(Plane(normal, normal * start_location));
  }
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_TPP_
