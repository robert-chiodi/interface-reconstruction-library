// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_SIMPLEX_WRAPPER_H_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_SIMPLEX_WRAPPER_H_

#include <array>

#include "src/generic_cutting/recursive_simplex_cutting/lookup_tables.h"
#include "src/geometry/general/geometry_type_traits.h"
#include "src/geometry/polygons/tri.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/helpers/SFINAE_boiler_plate.h"
#include "src/helpers/geometric_cutting_helpers.h"

namespace IRL {

template <class SimplexType, class Enable = void>
struct SimplexWrapper;

template <class VertexType>
class CutTetVertices;

template <class SimplexType>
struct SimplexWrapper<SimplexType, enable_if_t<is_tet<SimplexType>::value>> {
  static const constexpr bool isTet = true;
  static const constexpr bool isTri = false;

  static const constexpr LookupIndex_t simplex_nvert = 4;
  static const constexpr LookupIndex_t max_cut_simplex_nvert = 8;
  static constexpr LookupIndex_t
  numberOfSimplicesInVolumeBelowPlaneAfterCutting(
      const LookupIndex_t a_cutting_case) {
    return cut_tet_by_plane::number_of_negative_tets_after_cut[a_cutting_case];
  }
  static constexpr LookupIndex_t numberOfSimplicesInVolumeAfterCutting(
      const LookupIndex_t a_cutting_case) {
    return cut_tet_by_plane::number_of_tets_after_cut[a_cutting_case];
  }
  static double minimumAmountToTrack(void) {
    return global_constants::MINIMUM_VOLUME_TO_TRACK;
  }

  static constexpr bool isSimplexFullyBelowPlane(
      const LookupIndex_t a_cutting_case) {
    return a_cutting_case == 0;
  }

  static constexpr bool isSimplexFullyAbovePlane(
      const LookupIndex_t a_cutting_case) {
    return a_cutting_case == 15;
  }

  // Given a tet with signed distances to points,
  // returns the tet with vertices at intersection points
  /// \brief Calculates the vertices for a tet intersected with a plane.
  ///
  /// This function takes `a_tet` and the signed distance to each tet vertex,
  /// along with the unique integer identifier for the cutting configuration
  /// `a_cutting_case` to fill in a supplied Pt[8] array `a_vertices`. The
  /// first 4 elements in `a_vertices` are the original tet.
  /// The next vertices are vertices that lay on the intersection points
  /// between tet edges and the intersecting plane.
  /// These vertices are found using the function
  /// `ptPlaneIntersectsLine(...)`. The total number of vertices will be
  /// 4+`number_of_new_vertices_after_cut[a_cutting_case]`.
  ///
  /// \param[in] a_tet Tet that is being cut by plane and subdivided
  /// \param[in] a_vert_distances Double[4] array with the signed distance
  /// of each tet vertex from the cutting plane.
  /// \param[in] a_cutting_case Unique integer ID indicating
  /// cutting configuration.
  /// \param[out] a_vertices Array of up to 8 points which contains original
  /// tet vertices [0:3] and vertices at the intersection between tet edges
  /// and the cutting plane.
  static inline auto getSimplexVerticesAndIntersectionPoints(
      const SimplexType& a_tet,
      const std::array<double, simplex_nvert>& a_vert_distances,
      const LookupIndex_t a_cutting_case)
      -> CutTetVertices<typename SimplexType::pt_type>;

  /// \brief Returns a tet from tet vertices after intersection with a cutting
  /// plane.
  ///
  /// Given an array of `a_cut_tet_vertices` (where these vertices are the
  /// original tet vertices plus any vertices resulting from intersection with a
  /// cutting plane) the tet corresponding to `a_tet_number_to_get` in
  /// `verts_for_tets[][][]` from `lookup_tables.h` is returned.
  ///
  /// \param[in] a_cut_tet_vertices Array of vertices of a tet after
  /// intersection with a plane. Valid length of 1-8. \param[in] a_cutting_case
  /// The cutting case indicating signs of original tets vertices. \param[in]
  /// a_tet_number_to_get The tet to return from `verts_for_tets[][][]`.
  static inline constexpr auto simplexFromCutSimplexVertices(
      const CutTetVertices<typename SimplexType::pt_type>& a_cut_tet_vertices,
      const LookupIndex_t a_cutting_case,
      const LookupIndex_t a_tet_number_to_get)
      -> ProxyTet<CutTetVertices<typename SimplexType::pt_type>>;
};

template <class VertexType>
class CutTetVertices {
 public:
  using pt_type = VertexType;

  CutTetVertices(const pt_type& a_pt_0, const pt_type& a_pt_1,
                 const pt_type& a_pt_2, const pt_type& a_pt_3)
      : vertices_m{a_pt_0, a_pt_1, a_pt_2, a_pt_3} {}

  pt_type& operator[](const UnsignedIndex_t a_index) {
    assert(a_index < 8);
    return vertices_m[a_index];
  }

  const pt_type& operator[](const UnsignedIndex_t a_index) const {
    assert(a_index < 8);
    return vertices_m[a_index];
  }

 private:
  std::array<pt_type, 8> vertices_m;
};

template <class VertexType>
class CutTriangleVerticesAndPlane;

template <class SimplexType>
struct SimplexWrapper<SimplexType, enable_if_t<is_tri<SimplexType>::value>> {
  static const constexpr bool isTet = false;
  static const constexpr bool isTri = true;

  static const constexpr LookupIndex_t simplex_nvert = 3;
  static const constexpr LookupIndex_t max_cut_simplex_nvert = 5;
  static constexpr LookupIndex_t
  numberOfSimplicesInVolumeBelowPlaneAfterCutting(
      const LookupIndex_t a_cutting_case) {
    return cut_tri_by_plane::number_of_negative_tris_after_cut[a_cutting_case];
  }
  static constexpr UnsignedIndex_t numberOfSimplicesInVolumeAfterCutting(
      const LookupIndex_t a_cutting_case) {
    return cut_tri_by_plane::number_of_tris_after_cut[a_cutting_case];
  }
  static double minimumAmountToTrack(void) {
    return global_constants::MINIMUM_SURFACE_AREA_TO_TRACK;
  }

  static constexpr bool isSimplexFullyBelowPlane(
      const LookupIndex_t a_cutting_case) {
    return a_cutting_case == 0;
  }

  static constexpr bool isSimplexFullyAbovePlane(
      const LookupIndex_t a_cutting_case) {
    return a_cutting_case == 7;
  }

  static inline auto getSimplexVerticesAndIntersectionPoints(
      const SimplexType& a_tri,
      const std::array<double, simplex_nvert>& a_vert_distances,
      const LookupIndex_t a_cutting_case)
      -> CutTriangleVerticesAndPlane<typename SimplexType::pt_type>;

  static inline auto simplexFromCutSimplexVertices(
      const CutTriangleVerticesAndPlane<typename SimplexType::pt_type>&
          a_cut_tri_vertices,
      const LookupIndex_t a_cutting_case,
      const LookupIndex_t a_tri_number_to_get)
      -> ProxyTri<CutTriangleVerticesAndPlane<typename SimplexType::pt_type>>;
};

template <class VertexType>
class CutTriangleVerticesAndPlane {
 public:
  using pt_type = VertexType;

  CutTriangleVerticesAndPlane(void) : plane_of_existence_m(nullptr) {}

  explicit CutTriangleVerticesAndPlane(const pt_type& a_pt_0,
                                       const pt_type& a_pt_1,
                                       const pt_type& a_pt_2,
                                       const Plane* a_plane_of_existence)
      : vertices_m{a_pt_0, a_pt_1, a_pt_2},
        plane_of_existence_m(a_plane_of_existence) {}

  pt_type& operator[](const UnsignedIndex_t a_index) {
    return vertices_m[a_index];
  }

  const pt_type& operator[](const UnsignedIndex_t a_index) const {
    return vertices_m[a_index];
  }

  const Plane& getPlaneOfExistence(void) const { return *plane_of_existence_m; }

  ~CutTriangleVerticesAndPlane(void) = default;

 private:
  std::array<pt_type, 5> vertices_m;
  const Plane* plane_of_existence_m;
};

}  // namespace IRL

#include "src/generic_cutting/recursive_simplex_cutting/simplex_wrapper.tpp"

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_SIMPLEX_WRAPPER_H_
