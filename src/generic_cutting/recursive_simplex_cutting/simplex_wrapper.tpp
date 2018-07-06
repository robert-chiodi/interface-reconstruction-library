// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_SIMPLEX_WRAPPER_TPP_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_SIMPLEX_WRAPPER_TPP_

namespace IRL {

//******************************************************************* //
//     Function template definitions placed below this.
//******************************************************************* //
template <class SimplexType>
auto SimplexWrapper<SimplexType, enable_if_t<is_tet<SimplexType>::value>>::
    getSimplexVerticesAndIntersectionPoints(
        const SimplexType& a_tet,
        const std::array<double, simplex_nvert>& a_vert_distances,
        const LookupIndex_t a_cutting_case)
        -> CutTetVertices<typename SimplexType::pt_type> {
  assert(a_cutting_case < 16);
  using PtType = typename SimplexType::pt_type;

  CutTetVertices<typename SimplexType::pt_type>
      simplex_vertices_and_intersection_points{a_tet[0], a_tet[1], a_tet[2],
                                               a_tet[3]};
  for (LookupIndex_t v = 0;
       v < cut_tet_by_plane::number_of_new_vertices_after_cut[a_cutting_case];
       ++v) {
    auto v1 = cut_tet_by_plane::cut_vertices[a_cutting_case][0][v];
    auto v2 = cut_tet_by_plane::cut_vertices[a_cutting_case][1][v];
    assert(v1 < 4);
    assert(v2 < 4);
    assert(v1 != v2);
    simplex_vertices_and_intersection_points[static_cast<UnsignedIndex_t>(4 +
                                                                          v)] =
        a_tet[v1].fromEdgeIntersection(a_tet[v1], a_vert_distances[v1],
                                       a_tet[v2], a_vert_distances[v2]);
  }
  return simplex_vertices_and_intersection_points;
}

template <class SimplexType>
constexpr auto
SimplexWrapper<SimplexType, enable_if_t<is_tet<SimplexType>::value>>::
    simplexFromCutSimplexVertices(
        const CutTetVertices<typename SimplexType::pt_type>& a_cut_tet_vertices,
        const LookupIndex_t a_cutting_case,
        const LookupIndex_t a_tet_number_to_get)
        -> ProxyTet<CutTetVertices<typename SimplexType::pt_type>> {
  return {
      a_cut_tet_vertices,
      {cut_tet_by_plane::verts_for_tets[a_cutting_case][a_tet_number_to_get][0],
       cut_tet_by_plane::verts_for_tets[a_cutting_case][a_tet_number_to_get][1],
       cut_tet_by_plane::verts_for_tets[a_cutting_case][a_tet_number_to_get][2],
       cut_tet_by_plane::verts_for_tets[a_cutting_case][a_tet_number_to_get]
                                       [3]}};
}

template <class SimplexType>
inline auto
SimplexWrapper<SimplexType, enable_if_t<is_tri<SimplexType>::value>>::
    getSimplexVerticesAndIntersectionPoints(
        const SimplexType& a_tri,
        const std::array<double, simplex_nvert>& a_vert_distances,
        const LookupIndex_t a_cutting_case)
        -> CutTriangleVerticesAndPlane<typename SimplexType::pt_type> {
  assert(a_cutting_case < 8);
  using PtType = typename SimplexType::pt_type;

  CutTriangleVerticesAndPlane<typename SimplexType::pt_type> cut_triangle(
      a_tri[0], a_tri[1], a_tri[2], &a_tri.getPlaneOfExistence());
  for (LookupIndex_t v = 0;
       v < cut_tri_by_plane::number_of_new_vertices_after_cut[a_cutting_case];
       ++v) {
    auto v1 = cut_tri_by_plane::cut_vertices[a_cutting_case][0][v];
    auto v2 = cut_tri_by_plane::cut_vertices[a_cutting_case][1][v];
    assert(v1 < 3);
    assert(v2 < 3);
    cut_triangle[static_cast<UnsignedIndex_t>(3 + v)] =
        a_tri[v1].fromEdgeIntersection(a_tri[v1], a_vert_distances[v1],
                                       a_tri[v2], a_vert_distances[v2]);
  }
  return cut_triangle;
}

template <class SimplexType>
inline auto
SimplexWrapper<SimplexType, enable_if_t<is_tri<SimplexType>::value>>::
    simplexFromCutSimplexVertices(
        const CutTriangleVerticesAndPlane<typename SimplexType::pt_type>&
            a_cut_tri_vertices,
        const LookupIndex_t a_cutting_case,
        const LookupIndex_t a_tri_number_to_get)
        -> ProxyTri<
            CutTriangleVerticesAndPlane<typename SimplexType::pt_type>> {
  return {
      a_cut_tri_vertices,
      {cut_tri_by_plane::verts_for_tris[a_cutting_case][a_tri_number_to_get][0],
       cut_tri_by_plane::verts_for_tris[a_cutting_case][a_tri_number_to_get][1],
       cut_tri_by_plane::verts_for_tris[a_cutting_case][a_tri_number_to_get]
                                       [2]}};
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_SIMPLEX_WRAPPER_TPP_
