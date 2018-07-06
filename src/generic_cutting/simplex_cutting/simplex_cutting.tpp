// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_TPP_
#define SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_TPP_

#include <utility>

#include "src/generic_cutting/recursive_simplex_cutting/simplex_wrapper.h"

namespace IRL {
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
splitDecomposedPolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                        SegmentedHalfEdgePolyhedronType* a_clipped_polytope,
                        HalfEdgePolytopeType* a_complete_polytope,
                        const Plane& a_plane) {
  using SimplexType = typename SegmentedHalfEdgePolyhedronType::ProxyType;
  using pt_type = typename SegmentedHalfEdgePolyhedronType::pt_type;

  auto& vertex_list = a_polytope->getVertexList();
  auto& clipped_vertex_list = a_clipped_polytope->getVertexList();
  const auto starting_number_of_vertices = vertex_list.size();

  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    vertex_list.setDistance(n, a_plane.signedDistanceToPoint(vertex_list[n]));
  }

  UnsignedIndex_t full_gas_counter = 0;
  UnsignedIndex_t full_liquid_counter = 0;
  for (full_gas_counter = 0; full_gas_counter < starting_number_of_vertices;
       ++full_gas_counter) {
    if (vertex_list.getDistance(full_gas_counter) <= 0.0) {
      break;
    }
  }

  if (full_gas_counter == starting_number_of_vertices) {
    // Everything was above, nothing underneath
    std::swap(*a_clipped_polytope, *a_polytope);
    return;
  } else if (full_gas_counter == 0) {
    // Need to check if all liquid
    for (full_liquid_counter = 0;
         full_liquid_counter < starting_number_of_vertices;
         ++full_liquid_counter) {
      if (vertex_list.getDistance(full_liquid_counter) > 0.0) {
        break;
      }
    }
    if (full_liquid_counter == starting_number_of_vertices) {
      // Everything below, nothing above
      return;
    }
  }

  const auto original_last_tet_index =
      a_polytope->getNumberOfSimplicesInDecomposition() - 1;
  for (UnsignedIndex_t t = original_last_tet_index;
       t != static_cast<UnsignedIndex_t>(-1); --t) {
    const auto& current_simplex = a_polytope->getSimplexFromDecomposition(t);
    const auto& simplex_index_mapping = current_simplex.getIndexMapping();

    unsigned int cutting_case_uns = 0;
    for (UnsignedIndex_t n = 0; n < simplex_index_mapping.size(); ++n) {
      if (a_complete_polytope->getDistance(simplex_index_mapping[n]) > 0.0) {
        cutting_case_uns |= 1U << n;
      }
    }
    const LookupIndex_t cutting_case =
        static_cast<LookupIndex_t>(cutting_case_uns);

    if (SimplexWrapper<SimplexType>::isSimplexFullyBelowPlane(cutting_case)) {
      // Purely below plane, keep
      continue;
    } else if (SimplexWrapper<SimplexType>::isSimplexFullyAbovePlane(
                   cutting_case)) {
      // Purely clipped off, hand off to clipped one
      a_polytope->markAsNotValid(t);
      a_clipped_polytope->push_back(current_simplex);
      continue;
    }

    // If made it this far, an intersection has occurred.
    static std::array<UnsignedIndex_t,
                      SimplexWrapper<SimplexType>::max_cut_simplex_nvert>
        cut_simplex_vertex_indices;
    std::copy(simplex_index_mapping.begin(), simplex_index_mapping.end(),
              cut_simplex_vertex_indices.begin());
    const UnsignedIndex_t first_new_vertex_index = a_complete_polytope->size();
    a_complete_polytope->resize(
        first_new_vertex_index +
        cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (LookupIndex_t v = 0;
         v < cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case];
         ++v) {
      const auto v1 = cut_tet_by_plane::cut_vertices[cutting_case][0][v];
      const auto v2 = cut_tet_by_plane::cut_vertices[cutting_case][1][v];
      assert(v1 < 4);
      assert(v2 < 4);
      assert(v1 != v2);
      (*a_complete_polytope)[first_new_vertex_index + v] =
          current_simplex[v1].fromEdgeIntersection(
              current_simplex[v1],
              a_complete_polytope->getDistance(simplex_index_mapping[v1]),
              current_simplex[v2],
              a_complete_polytope->getDistance(simplex_index_mapping[v2]));
    }

    const UnsignedIndex_t current_vertex_list_size = vertex_list.size();
    vertex_list.resize(
        current_vertex_list_size +
        cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    const UnsignedIndex_t current_clipped_vertex_list_size =
        clipped_vertex_list.size();
    clipped_vertex_list.resize(
        current_clipped_vertex_list_size +
        cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (UnsignedIndex_t v = 0;
         v <
         static_cast<UnsignedIndex_t>(
             cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case]);
         ++v) {
      cut_simplex_vertex_indices[4 + v] = first_new_vertex_index + v;
      vertex_list.getVertexIndex(current_vertex_list_size + v) =
          first_new_vertex_index + v;
      clipped_vertex_list.getVertexIndex(current_clipped_vertex_list_size + v) =
          first_new_vertex_index + v;
    }

    // Create and add new proxy simplices underneath plane
    for (LookupIndex_t tet_below = 0;
         tet_below <
         SimplexWrapper<SimplexType>::
             numberOfSimplicesInVolumeBelowPlaneAfterCutting(cutting_case);
         ++tet_below) {
      const auto simplex_from_volume_below_plane = SimplexType(
          *a_complete_polytope,
          {cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][0]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][1]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][2]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][3]]});

      if (simplex_from_volume_below_plane.calculateAbsoluteVolume() >
          SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
        a_polytope->push_back(simplex_from_volume_below_plane);
      }
    }

    // Create and add new above plane simplices
    for (LookupIndex_t tet_above = SimplexWrapper<SimplexType>::
             numberOfSimplicesInVolumeBelowPlaneAfterCutting(cutting_case);
         tet_above <
         SimplexWrapper<SimplexType>::numberOfSimplicesInVolumeAfterCutting(
             cutting_case);
         ++tet_above) {
      const auto simplex_from_volume_above_plane = SimplexType(
          *a_complete_polytope,
          {cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_above][0]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_above][1]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_above][2]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_above][3]]});

      if (simplex_from_volume_above_plane.calculateAbsoluteVolume() >
          SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
        a_clipped_polytope->push_back(simplex_from_volume_above_plane);
      }
    }

    // Remove the simplex at t
    a_polytope->markAsNotValid(t);
  }

  // Clean up
  UnsignedIndex_t number_kept = 0;
  for (UnsignedIndex_t n = 0; n < original_last_tet_index + 1; ++n) {
    if (a_polytope->isValid(n)) {
      a_polytope->getSimplexFromDecomposition(number_kept) =
          a_polytope->getSimplexFromDecomposition(n);
      ++number_kept;
    }
  }
  for (UnsignedIndex_t n = original_last_tet_index + 1;
       n < a_polytope->getNumberOfSimplicesInDecomposition(); ++n) {
    a_polytope->getSimplexFromDecomposition(
        number_kept + (n - (original_last_tet_index + 1))) =
        a_polytope->getSimplexFromDecomposition(n);
  }
  number_kept = number_kept +
                a_polytope->getNumberOfSimplicesInDecomposition() -
                (original_last_tet_index + 1);
  a_polytope->setNumberOfSimplicesInDecomposition(number_kept);
  a_polytope->setNumberOfValidSimplices(number_kept);
  for (UnsignedIndex_t n = 0; n < number_kept; ++n) {
    a_polytope->markAsValid(n);
  }

  number_kept = 0;
  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    if (vertex_list.getDistance(n) > 0.0) {
      clipped_vertex_list.push_back(vertex_list.getVertexIndex(n));
    } else {
      vertex_list.getVertexIndex(number_kept) = vertex_list.getVertexIndex(n);
      ++number_kept;
    }
  }
  const auto new_length = vertex_list.size();
  for (UnsignedIndex_t n = starting_number_of_vertices; n < new_length; ++n) {
    vertex_list.getVertexIndex(number_kept +
                               (n - starting_number_of_vertices)) =
        vertex_list.getVertexIndex(n);
  }
  vertex_list.resize(number_kept + (new_length - starting_number_of_vertices));
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
// Splitting a simplex-segmented polygon by a plane
////////////////////////////////////////////////////
////////////////////////////////////////////////////
// During generation of a SegmentedDecomposedPolygon, the plane that the
// polygon exists on should have been created. Any polygon split from this
// polygon should also exist on this same plane, since that's how polygons
// work (on a 2D plane). This plane is used during the calculation of the
// moments in order to assign a consistent reference direction around which
// to determine signed normals.
template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
splitDecomposedPolytope(SegmentedHalfEdgePolygonType* a_polytope,
                        SegmentedHalfEdgePolygonType* a_clipped_polytope,
                        HalfEdgePolytopeType* a_complete_polytope,
                        const Plane& a_plane) {
  using SimplexType = typename SegmentedHalfEdgePolygonType::ProxyType;
  using pt_type = typename SegmentedHalfEdgePolygonType::pt_type;

  a_clipped_polytope->setPlaneOfExistence(&a_polytope->getPlaneOfExistence());

  auto& vertex_list = a_polytope->getVertexList();
  auto& clipped_vertex_list = a_clipped_polytope->getVertexList();
  const auto starting_number_of_vertices = vertex_list.size();

  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    vertex_list.setDistance(n, a_plane.signedDistanceToPoint(vertex_list[n]));
  }

  UnsignedIndex_t full_gas_counter = 0;
  UnsignedIndex_t full_liquid_counter = 0;
  for (full_gas_counter = 0; full_gas_counter < starting_number_of_vertices;
       ++full_gas_counter) {
    if (vertex_list.getDistance(full_gas_counter) <= 0.0) {
      break;
    }
  }

  if (full_gas_counter == starting_number_of_vertices) {
    // Everything was above, nothing underneath
    std::swap(*a_clipped_polytope, *a_polytope);
    return;
  } else if (full_gas_counter == 0) {
    // Need to check if all liquid
    for (full_liquid_counter = 0;
         full_liquid_counter < starting_number_of_vertices;
         ++full_liquid_counter) {
      if (vertex_list.getDistance(full_liquid_counter) > 0.0) {
        break;
      }
    }
    if (full_liquid_counter == starting_number_of_vertices) {
      // Everything below, nothing above
      return;
    }
  }

  const auto original_last_tet_index =
      a_polytope->getNumberOfSimplicesInDecomposition() - 1;
  for (UnsignedIndex_t t = original_last_tet_index;
       t != static_cast<UnsignedIndex_t>(-1); --t) {
    const auto& current_simplex = a_polytope->getSimplexFromDecomposition(t);
    const auto& simplex_index_mapping = current_simplex.getIndexMapping();

    unsigned int cutting_case_uns = 0;
    for (UnsignedIndex_t n = 0; n < simplex_index_mapping.size(); ++n) {
      if (a_complete_polytope->getDistance(simplex_index_mapping[n]) > 0.0) {
        cutting_case_uns |= 1U << n;
      }
    }
    const LookupIndex_t cutting_case =
        static_cast<LookupIndex_t>(cutting_case_uns);

    if (SimplexWrapper<SimplexType>::isSimplexFullyBelowPlane(cutting_case)) {
      // Purely below plane, keep
      continue;
    } else if (SimplexWrapper<SimplexType>::isSimplexFullyAbovePlane(
                   cutting_case)) {
      // Purely clipped off, hand off to clipped one
      a_polytope->markAsNotValid(t);
      a_clipped_polytope->push_back(current_simplex);
      continue;
    }

    // If made it this far, an intersection has occurred.
    static std::array<UnsignedIndex_t,
                      SimplexWrapper<SimplexType>::max_cut_simplex_nvert>
        cut_simplex_vertex_indices;
    std::copy(simplex_index_mapping.begin(), simplex_index_mapping.end(),
              cut_simplex_vertex_indices.begin());
    const UnsignedIndex_t first_new_vertex_index = a_complete_polytope->size();
    a_complete_polytope->resize(
        first_new_vertex_index +
        cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (LookupIndex_t v = 0;
         v < cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case];
         ++v) {
      const auto v1 = cut_tri_by_plane::cut_vertices[cutting_case][0][v];
      const auto v2 = cut_tri_by_plane::cut_vertices[cutting_case][1][v];
      assert(v1 < 3);
      assert(v2 < 3);
      assert(v1 != v2);
      (*a_complete_polytope)[first_new_vertex_index + v] =
          current_simplex[v1].fromEdgeIntersection(
              current_simplex[v1],
              a_complete_polytope->getDistance(simplex_index_mapping[v1]),
              current_simplex[v2],
              a_complete_polytope->getDistance(simplex_index_mapping[v2]));
    }

    const UnsignedIndex_t current_vertex_list_size = vertex_list.size();
    vertex_list.resize(
        current_vertex_list_size +
        cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    const UnsignedIndex_t current_clipped_vertex_list_size =
        clipped_vertex_list.size();
    clipped_vertex_list.resize(
        current_clipped_vertex_list_size +
        cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (UnsignedIndex_t v = 0;
         v <
         static_cast<UnsignedIndex_t>(
             cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case]);
         ++v) {
      cut_simplex_vertex_indices[3 + v] = first_new_vertex_index + v;
      vertex_list.getVertexIndex(current_vertex_list_size + v) =
          first_new_vertex_index + v;
      clipped_vertex_list.getVertexIndex(current_clipped_vertex_list_size + v) =
          first_new_vertex_index + v;
    }

    // Create and add new proxy simplices underneath plane
    for (LookupIndex_t tet_below = 0;
         tet_below <
         SimplexWrapper<SimplexType>::
             numberOfSimplicesInVolumeBelowPlaneAfterCutting(cutting_case);
         ++tet_below) {
      const auto simplex_from_volume_below_plane = SimplexType(
          *a_complete_polytope,
          {cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_below][0]],
           cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_below][1]],
           cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_below][2]]});

      if (simplex_from_volume_below_plane.calculateAbsoluteVolume() >
          SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
        a_polytope->push_back(simplex_from_volume_below_plane);
      }
    }

    // Create and add new above plane simplices
    for (LookupIndex_t tet_above = SimplexWrapper<SimplexType>::
             numberOfSimplicesInVolumeBelowPlaneAfterCutting(cutting_case);
         tet_above <
         SimplexWrapper<SimplexType>::numberOfSimplicesInVolumeAfterCutting(
             cutting_case);
         ++tet_above) {
      const auto simplex_from_volume_above_plane = SimplexType(
          *a_complete_polytope,
          {cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_above][0]],
           cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_above][1]],
           cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_above][2]]});

      if (simplex_from_volume_above_plane.calculateAbsoluteVolume() >
          SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
        a_clipped_polytope->push_back(simplex_from_volume_above_plane);
      }
    }

    // Remove the simplex at t
    a_polytope->markAsNotValid(t);
  }

  // Clean up
  UnsignedIndex_t number_kept = 0;
  for (UnsignedIndex_t n = 0; n < original_last_tet_index + 1; ++n) {
    if (a_polytope->isValid(n)) {
      a_polytope->getSimplexFromDecomposition(number_kept) =
          a_polytope->getSimplexFromDecomposition(n);
      ++number_kept;
    }
  }
  for (UnsignedIndex_t n = original_last_tet_index + 1;
       n < a_polytope->getNumberOfSimplicesInDecomposition(); ++n) {
    a_polytope->getSimplexFromDecomposition(
        number_kept + (n - (original_last_tet_index + 1))) =
        a_polytope->getSimplexFromDecomposition(n);
  }
  number_kept = number_kept +
                a_polytope->getNumberOfSimplicesInDecomposition() -
                (original_last_tet_index + 1);
  a_polytope->setNumberOfSimplicesInDecomposition(number_kept);
  a_polytope->setNumberOfValidSimplices(number_kept);
  for (UnsignedIndex_t n = 0; n < number_kept; ++n) {
    a_polytope->markAsValid(n);
  }

  number_kept = 0;
  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    if (vertex_list.getDistance(n) > 0.0) {
      clipped_vertex_list.push_back(vertex_list.getVertexIndex(n));
    } else {
      vertex_list.getVertexIndex(number_kept) = vertex_list.getVertexIndex(n);
      ++number_kept;
    }
  }
  const auto new_length = vertex_list.size();
  for (UnsignedIndex_t n = starting_number_of_vertices; n < new_length; ++n) {
    vertex_list.getVertexIndex(number_kept +
                               (n - starting_number_of_vertices)) =
        vertex_list.getVertexIndex(n);
  }
  vertex_list.resize(number_kept + (new_length - starting_number_of_vertices));
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Truncating instead of splitting below this ///
//////////////////////////////////////////////////
//////////////////////////////////////////////////
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
truncateDecomposedPolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                           HalfEdgePolytopeType* a_complete_polytope,
                           const Plane& a_plane) {
  using SimplexType = typename SegmentedHalfEdgePolyhedronType::ProxyType;
  using pt_type = typename SegmentedHalfEdgePolyhedronType::pt_type;

  auto& vertex_list = a_polytope->getVertexList();
  const auto starting_number_of_vertices = vertex_list.size();

  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    vertex_list.setDistance(n, a_plane.signedDistanceToPoint(vertex_list[n]));
  }

  UnsignedIndex_t full_gas_counter = 0;
  UnsignedIndex_t full_liquid_counter = 0;
  for (full_gas_counter = 0; full_gas_counter < starting_number_of_vertices;
       ++full_gas_counter) {
    if (vertex_list.getDistance(full_gas_counter) <= 0.0) {
      break;
    }
  }

  if (full_gas_counter == starting_number_of_vertices) {
    // Everything was above, nothing underneath
    // Everything was above, nothing underneath
    a_polytope->setNumberOfSimplicesInDecomposition(0);
    a_polytope->setNumberOfValidSimplices(0);
    vertex_list.resize(0);
    return;
    return;
  } else if (full_gas_counter == 0) {
    // Need to check if all liquid
    for (full_liquid_counter = 0;
         full_liquid_counter < starting_number_of_vertices;
         ++full_liquid_counter) {
      if (vertex_list.getDistance(full_liquid_counter) > 0.0) {
        break;
      }
    }
    if (full_liquid_counter == starting_number_of_vertices) {
      // Everything below, nothing above
      return;
    }
  }

  const auto original_last_tet_index =
      a_polytope->getNumberOfSimplicesInDecomposition() - 1;
  for (UnsignedIndex_t t = original_last_tet_index;
       t != static_cast<UnsignedIndex_t>(-1); --t) {
    const auto& current_simplex = a_polytope->getSimplexFromDecomposition(t);
    const auto& simplex_index_mapping = current_simplex.getIndexMapping();

    unsigned int cutting_case_uns = 0;
    for (UnsignedIndex_t n = 0; n < simplex_index_mapping.size(); ++n) {
      if (a_complete_polytope->getDistance(simplex_index_mapping[n]) > 0.0) {
        cutting_case_uns |= 1U << n;
      }
    }
    const LookupIndex_t cutting_case =
        static_cast<LookupIndex_t>(cutting_case_uns);

    if (SimplexWrapper<SimplexType>::isSimplexFullyBelowPlane(cutting_case)) {
      // Purely below plane, keep
      continue;
    } else if (SimplexWrapper<SimplexType>::isSimplexFullyAbovePlane(
                   cutting_case)) {
      // Purely clipped off, mark as invalid
      a_polytope->markAsNotValid(t);
      continue;
    }

    // If made it this far, an intersection has occurred.
    static std::array<UnsignedIndex_t,
                      SimplexWrapper<SimplexType>::max_cut_simplex_nvert>
        cut_simplex_vertex_indices;
    std::copy(simplex_index_mapping.begin(), simplex_index_mapping.end(),
              cut_simplex_vertex_indices.begin());
    const UnsignedIndex_t first_new_vertex_index = a_complete_polytope->size();
    a_complete_polytope->resize(
        first_new_vertex_index +
        cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (LookupIndex_t v = 0;
         v < cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case];
         ++v) {
      const auto v1 = cut_tet_by_plane::cut_vertices[cutting_case][0][v];
      const auto v2 = cut_tet_by_plane::cut_vertices[cutting_case][1][v];
      assert(v1 < 4);
      assert(v2 < 4);
      assert(v1 != v2);
      (*a_complete_polytope)[first_new_vertex_index + v] =
          current_simplex[v1].fromEdgeIntersection(
              current_simplex[v1],
              a_complete_polytope->getDistance(simplex_index_mapping[v1]),
              current_simplex[v2],
              a_complete_polytope->getDistance(simplex_index_mapping[v2]));
    }

    const UnsignedIndex_t current_vertex_list_size = vertex_list.size();
    vertex_list.resize(
        current_vertex_list_size +
        cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (UnsignedIndex_t v = 0;
         v <
         static_cast<UnsignedIndex_t>(
             cut_tet_by_plane::number_of_new_vertices_after_cut[cutting_case]);
         ++v) {
      cut_simplex_vertex_indices[4 + v] = first_new_vertex_index + v;
      vertex_list.getVertexIndex(current_vertex_list_size + v) =
          first_new_vertex_index + v;
    }

    // Create and add new proxy simplices underneath plane
    for (LookupIndex_t tet_below = 0;
         tet_below <
         SimplexWrapper<SimplexType>::
             numberOfSimplicesInVolumeBelowPlaneAfterCutting(cutting_case);
         ++tet_below) {
      const auto simplex_from_volume_below_plane = SimplexType(
          *a_complete_polytope,
          {cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][0]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][1]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][2]],
           cut_simplex_vertex_indices
               [cut_tet_by_plane::verts_for_tets[cutting_case][tet_below][3]]});

      if (simplex_from_volume_below_plane.calculateAbsoluteVolume() >
          SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
        a_polytope->push_back(simplex_from_volume_below_plane);
      }
    }

    // Remove the simplex at t
    a_polytope->markAsNotValid(t);
  }

  // Clean up the stored Simplices
  UnsignedIndex_t number_kept = 0;
  for (UnsignedIndex_t n = 0; n < original_last_tet_index + 1; ++n) {
    if (a_polytope->isValid(n)) {
      a_polytope->getSimplexFromDecomposition(number_kept) =
          a_polytope->getSimplexFromDecomposition(n);
      ++number_kept;
    }
  }
  for (UnsignedIndex_t n = original_last_tet_index + 1;
       n < a_polytope->getNumberOfSimplicesInDecomposition(); ++n) {
    a_polytope->getSimplexFromDecomposition(
        number_kept + (n - (original_last_tet_index + 1))) =
        a_polytope->getSimplexFromDecomposition(n);
  }
  number_kept = number_kept +
                a_polytope->getNumberOfSimplicesInDecomposition() -
                (original_last_tet_index + 1);
  a_polytope->setNumberOfSimplicesInDecomposition(number_kept);
  a_polytope->setNumberOfValidSimplices(number_kept);
  for (UnsignedIndex_t n = 0; n < number_kept; ++n) {
    a_polytope->markAsValid(n);
  }

  // Clean up the stored vertices
  number_kept = 0;
  const UnsignedIndex_t new_length = vertex_list.size();
  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    if (vertex_list.getDistance(n) <= 0.0) {
      vertex_list.getVertexIndex(number_kept) = vertex_list.getVertexIndex(n);
      ++number_kept;
    }
  }
  for (UnsignedIndex_t n = starting_number_of_vertices; n < new_length; ++n) {
    vertex_list.getVertexIndex(number_kept +
                               (n - starting_number_of_vertices)) =
        vertex_list.getVertexIndex(n);
  }
  vertex_list.resize(number_kept + (new_length - starting_number_of_vertices));
}

////////////////////////////////////////
// Triangle truncation below this     //
////////////////////////////////////////
template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
truncateDecomposedPolytope(SegmentedHalfEdgePolygonType* a_polytope,
                           HalfEdgePolytopeType* a_complete_polytope,
                           const Plane& a_plane) {
  using SimplexType = typename SegmentedHalfEdgePolygonType::ProxyType;
  using pt_type = typename SegmentedHalfEdgePolygonType::pt_type;

  auto& vertex_list = a_polytope->getVertexList();
  const auto starting_number_of_vertices = vertex_list.size();

  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    vertex_list.setDistance(n, a_plane.signedDistanceToPoint(vertex_list[n]));
  }

  UnsignedIndex_t full_gas_counter = 0;
  UnsignedIndex_t full_liquid_counter = 0;
  for (full_gas_counter = 0; full_gas_counter < starting_number_of_vertices;
       ++full_gas_counter) {
    if (vertex_list.getDistance(full_gas_counter) <= 0.0) {
      break;
    }
  }

  if (full_gas_counter == starting_number_of_vertices) {
    // Everything was above, nothing underneath
    // Everything was above, nothing underneath
    a_polytope->setNumberOfSimplicesInDecomposition(0);
    a_polytope->setNumberOfValidSimplices(0);
    vertex_list.resize(0);
    return;
    return;
  } else if (full_gas_counter == 0) {
    // Need to check if all liquid
    for (full_liquid_counter = 0;
         full_liquid_counter < starting_number_of_vertices;
         ++full_liquid_counter) {
      if (vertex_list.getDistance(full_liquid_counter) > 0.0) {
        break;
      }
    }
    if (full_liquid_counter == starting_number_of_vertices) {
      // Everything below, nothing above
      return;
    }
  }

  const auto original_last_tet_index =
      a_polytope->getNumberOfSimplicesInDecomposition() - 1;
  for (UnsignedIndex_t t = original_last_tet_index;
       t != static_cast<UnsignedIndex_t>(-1); --t) {
    const auto& current_simplex = a_polytope->getSimplexFromDecomposition(t);
    const auto& simplex_index_mapping = current_simplex.getIndexMapping();

    unsigned int cutting_case_uns = 0;
    for (UnsignedIndex_t n = 0; n < simplex_index_mapping.size(); ++n) {
      if (a_complete_polytope->getDistance(simplex_index_mapping[n]) > 0.0) {
        cutting_case_uns |= 1U << n;
      }
    }
    const LookupIndex_t cutting_case =
        static_cast<LookupIndex_t>(cutting_case_uns);

    if (SimplexWrapper<SimplexType>::isSimplexFullyBelowPlane(cutting_case)) {
      // Purely below plane, keep
      continue;
    } else if (SimplexWrapper<SimplexType>::isSimplexFullyAbovePlane(
                   cutting_case)) {
      // Purely clipped off, hand off to clipped one
      a_polytope->markAsNotValid(t);
      continue;
    }

    // If made it this far, an intersection has occurred.
    static std::array<UnsignedIndex_t,
                      SimplexWrapper<SimplexType>::max_cut_simplex_nvert>
        cut_simplex_vertex_indices;
    std::copy(simplex_index_mapping.begin(), simplex_index_mapping.end(),
              cut_simplex_vertex_indices.begin());
    const UnsignedIndex_t first_new_vertex_index = a_complete_polytope->size();
    a_complete_polytope->resize(
        first_new_vertex_index +
        cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (LookupIndex_t v = 0;
         v < cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case];
         ++v) {
      const auto v1 = cut_tri_by_plane::cut_vertices[cutting_case][0][v];
      const auto v2 = cut_tri_by_plane::cut_vertices[cutting_case][1][v];
      assert(v1 < 3);
      assert(v2 < 3);
      assert(v1 != v2);
      (*a_complete_polytope)[first_new_vertex_index + v] =
          current_simplex[v1].fromEdgeIntersection(
              current_simplex[v1],
              a_complete_polytope->getDistance(simplex_index_mapping[v1]),
              current_simplex[v2],
              a_complete_polytope->getDistance(simplex_index_mapping[v2]));
    }

    const UnsignedIndex_t current_vertex_list_size = vertex_list.size();
    vertex_list.resize(
        current_vertex_list_size +
        cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case]);
    for (UnsignedIndex_t v = 0;
         v <
         static_cast<UnsignedIndex_t>(
             cut_tri_by_plane::number_of_new_vertices_after_cut[cutting_case]);
         ++v) {
      cut_simplex_vertex_indices[3 + v] = first_new_vertex_index + v;
      vertex_list.getVertexIndex(current_vertex_list_size + v) =
          first_new_vertex_index + v;
    }

    // Create and add new proxy simplices underneath plane
    for (LookupIndex_t tet_below = 0;
         tet_below <
         SimplexWrapper<SimplexType>::
             numberOfSimplicesInVolumeBelowPlaneAfterCutting(cutting_case);
         ++tet_below) {
      const auto simplex_from_volume_below_plane = SimplexType(
          *a_complete_polytope,
          {cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_below][0]],
           cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_below][1]],
           cut_simplex_vertex_indices
               [cut_tri_by_plane::verts_for_tris[cutting_case][tet_below][2]]});

      if (simplex_from_volume_below_plane.calculateAbsoluteVolume() >
          SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
        a_polytope->push_back(simplex_from_volume_below_plane);
      }
    }

    // Remove the simplex at t
    a_polytope->markAsNotValid(t);
  }

  // Clean up
  UnsignedIndex_t number_kept = 0;
  for (UnsignedIndex_t n = 0; n < original_last_tet_index + 1; ++n) {
    if (a_polytope->isValid(n)) {
      a_polytope->getSimplexFromDecomposition(number_kept) =
          a_polytope->getSimplexFromDecomposition(n);
      ++number_kept;
    }
  }
  for (UnsignedIndex_t n = original_last_tet_index + 1;
       n < a_polytope->getNumberOfSimplicesInDecomposition(); ++n) {
    a_polytope->getSimplexFromDecomposition(
        number_kept + (n - (original_last_tet_index + 1))) =
        a_polytope->getSimplexFromDecomposition(n);
  }
  number_kept = number_kept +
                a_polytope->getNumberOfSimplicesInDecomposition() -
                (original_last_tet_index + 1);
  a_polytope->setNumberOfSimplicesInDecomposition(number_kept);
  a_polytope->setNumberOfValidSimplices(number_kept);
  for (UnsignedIndex_t n = 0; n < number_kept; ++n) {
    a_polytope->markAsValid(n);
  }

  number_kept = 0;
  const UnsignedIndex_t new_length = vertex_list.size();
  for (UnsignedIndex_t n = 0; n < starting_number_of_vertices; ++n) {
    if (vertex_list.getDistance(n) <= 0.0) {
      vertex_list.getVertexIndex(number_kept) = vertex_list.getVertexIndex(n);
      ++number_kept;
    }
  }
  for (UnsignedIndex_t n = starting_number_of_vertices; n < new_length; ++n) {
    vertex_list.getVertexIndex(number_kept +
                               (n - starting_number_of_vertices)) =
        vertex_list.getVertexIndex(n);
  }
  vertex_list.resize(number_kept + (new_length - starting_number_of_vertices));
}
}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_TPP_
