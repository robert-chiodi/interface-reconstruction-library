// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_TPP_
#define SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_TPP_

namespace IRL {

// Foward declare getVolumeMoments to avoid circular dependency.
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
__attribute__((hot)) inline ReturnType getVolumeMoments(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction);

/// Setting the half edge structure using lazy evaluation.
/// During the first access by the particular object type,
/// the HalfEdgePolyhedron will be generated in
/// complete_polyhedron_buffer using the
/// classes method and saved to a separate
/// HalfEdgePolyhedron. Since the half edge structure
/// for each particular object will remain the same, and since
/// the memory for complete_polyhedron_buffer is static for the entire
/// execution of the program, the object can now be recreated
/// on all subsequent accesses through direct copying of the
/// original half edge structure stored in the copy made during the
/// first accessing of this function. Then only the vertex locations
/// need to be reset for each new instance of the class object. For a
/// polygon, the plane of existence will also need to be reset
/// during each access.

template <class VertexType>
HalfEdgePolyhedron<VertexType>& getHalfEdgePolyhedronStorage(void);

template <class VertexType>
HalfEdgePolygon<VertexType>& getHalfEdgePolygonStorage(void);

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry);

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry);

template <class EncompassingGeometryType>
enable_if_t<is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry);

template <class EncompassingGeometryType>
enable_if_t<is_polygon<EncompassingGeometryType>::value &&
                !is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry);

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolyhedron<typename PolytopeType::face_type,
                                        typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType* a_complete_polytope);

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolyhedron<typename PolytopeType::face_type,
                                        typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType* a_complete_polytope);

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polygon<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolygon<typename PolytopeType::face_type,
                                     typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType* a_complete_polytope);

//******************************************************************* //
//     Function template definitions placed below this
//******************************************************************* //
template <class ReturnType, class EncompassingType, class ReconstructionType>
ReturnType cutThroughHalfEdgeStructures(
    const EncompassingType& a_polytope,
    const ReconstructionType& a_reconstruction) {
  auto& complete_polytope = setHalfEdgeStructure(a_polytope);
  auto half_edge_polytope =
      generateSegmentedVersion<EncompassingType>(&complete_polytope);
  assert(half_edge_polytope.checkValidHalfEdgeStructure());
  return getVolumeMoments<ReturnType, HalfEdgeCutting>(
      &half_edge_polytope, &complete_polytope, a_reconstruction);
}

template <class VertexType>
HalfEdgePolyhedron<VertexType>& getHalfEdgePolyhedronStorage(void) {
  static HalfEdgePolyhedron<VertexType> storage;
  return storage;
}

template <class VertexType>
HalfEdgePolygon<VertexType>& getHalfEdgePolygonStorage(void) {
  static HalfEdgePolygon<VertexType> storage;
  return storage;
}

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry) {
  static bool already_set = false;
  static HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>
      half_edge_geometry_template;

  HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>&
      complete_polyhedron_buffer = getHalfEdgePolyhedronStorage<
          typename EncompassingGeometryType::pt_type>();

  if (already_set) {
    complete_polyhedron_buffer = half_edge_geometry_template;
    complete_polyhedron_buffer.setVertexLocations(a_geometry);
  } else {
    already_set = true;
    a_geometry.setHalfEdgeVersion(&complete_polyhedron_buffer);
    half_edge_geometry_template = complete_polyhedron_buffer;
  }
  return complete_polyhedron_buffer;
}

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry) {
  // Can't lazy evaluate GeneralPolyhedron's because they change, i.e.
  // since a GeneralPolyhedron doesn't have a fixed number of vertices or
  // connectivity like other polyhedrons in IRL, the connectivity cannot simply
  // be stored and copied over.
  HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>&
      complete_polyhedron_buffer = getHalfEdgePolyhedronStorage<
          typename EncompassingGeometryType::pt_type>();

  a_geometry.setHalfEdgeVersion(&complete_polyhedron_buffer);
  return complete_polyhedron_buffer;
}

template <class EncompassingGeometryType>
enable_if_t<is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry) {
  static bool already_set = false;
  static HalfEdgePolygon<typename EncompassingGeometryType::pt_type>
      half_edge_geometry_template;

  HalfEdgePolygon<
      typename EncompassingGeometryType::pt_type>& complete_polygon_buffer =
      getHalfEdgePolygonStorage<typename EncompassingGeometryType::pt_type>();

  if (already_set) {
    complete_polygon_buffer = half_edge_geometry_template;
  } else {
    already_set = true;
    a_geometry.setHalfEdgeVersion(&complete_polygon_buffer);
    half_edge_geometry_template = complete_polygon_buffer;
  }
  complete_polygon_buffer.setVertexLocations(a_geometry);
  complete_polygon_buffer.setPlaneOfExistence(a_geometry.getPlaneOfExistence());
  return complete_polygon_buffer;
}

template <class EncompassingGeometryType>
enable_if_t<is_polygon<EncompassingGeometryType>::value &&
                !is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type>&>
setHalfEdgeStructure(const EncompassingGeometryType& a_geometry) {
  // Can't lazy evaluate polygons because they change, i.e.
  // since a Polygon doesn't have a fixed number of vertices like
  // a polyhedron does in IRL, the connectivity cannot simply
  // be stored and copied over.
  HalfEdgePolygon<
      typename EncompassingGeometryType::pt_type>& complete_polygon_buffer =
      getHalfEdgePolygonStorage<typename EncompassingGeometryType::pt_type>();
  a_geometry.setHalfEdgeVersion(&complete_polygon_buffer);
  return complete_polygon_buffer;
}

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolyhedron<typename PolytopeType::face_type,
                                        typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType* a_complete_polytope) {
  static bool already_set = false;
  static SegmentedHalfEdgePolyhedron<typename PolytopeType::face_type,
                                     typename PolytopeType::vertex_type>
      segmented_half_edge_geometry_template;

  if (!already_set) {
    already_set = true;
    segmented_half_edge_geometry_template =
        a_complete_polytope->generateSegmentedPolyhedron();
  }
  return segmented_half_edge_geometry_template;
}

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolyhedron<typename PolytopeType::face_type,
                                        typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType* a_complete_polytope) {
  return a_complete_polytope->generateSegmentedPolyhedron();
}

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polygon<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolygon<typename PolytopeType::face_type,
                                     typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType* a_complete_polytope) {
  auto segmented_half_edge_geometry_template =
      a_complete_polytope->generateSegmentedPolygon();
  segmented_half_edge_geometry_template.setPlaneOfExistence(
      &a_complete_polytope->getPlaneOfExistence());
  return segmented_half_edge_geometry_template;
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_TPP_
