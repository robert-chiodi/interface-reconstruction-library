// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_TPP_
#define IRL_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_TPP_

namespace IRL {

// Foward declare getVolumeMoments to avoid circular dependency.
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
__attribute__((hot)) inline ReturnType
getVolumeMoments(SegmentedPolytopeType *a_polytope,
                 HalfEdgePolytopeType *a_complete_polytope,
                 const ReconstructionType &a_reconstruction);

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
template <class PolytopeType> class HalfEdgeGeometryInitializer {

public:
  HalfEdgeGeometryInitializer(void) = default;

  void setAllToUpdate(void) {
    std::fill(needs_updating_m.begin(), needs_updating_m.end(), 1);
  }

  void setAllOthersToUpdate(std::size_t a_type_id) {
    this->setAllToUpdate();
    assert(a_type_id < needs_updating_m.size());
    needs_updating_m[a_type_id] = 0;
  }

  template <class GeometryType>
  std::size_t getID(const GeometryType &a_geometry) {
    static bool already_set = false;
    static std::size_t type_id = static_cast<std::size_t>(-1);
    if (!already_set) {
      already_set = true;
      type_id = template_polytopes_m.size();
      const auto starting_capacity = central_polytope_storage_m.rawCapacity();
      PolytopeType tmp;
      template_polytopes_m.push_back(PolytopeType());
      needs_updating_m.push_back(0);
      this->resetTemplatePolytope(a_geometry, type_id);
      if (central_polytope_storage_m.rawCapacity() > starting_capacity) {
        // New type forced a reallocation. Invalidates all other
        // polytope templates.
        this->setAllOthersToUpdate(type_id);
      }
    }
    return type_id;
  }

  void resetCentralStorageToCurrentSize(void) {
    const auto current_capacity = central_polytope_storage_m.rawCapacity();
    if (central_polytope_storage_m.firstBlockCapacity() < current_capacity) {
      // There are additional blocks added on. Resize to have one contiguous
      // block
      central_polytope_storage_m.resizeRawCapacity(current_capacity);
      this->setAllToUpdate();
    }
  }

  template <class GeometryType>
  void resetTemplatePolytope(const GeometryType &a_geometry,
                             const std::size_t a_type_id) {
    assert(a_type_id < template_polytopes_m.size());
    assert(a_type_id < needs_updating_m.size());
    a_geometry.setHalfEdgeVersion(&central_polytope_storage_m);
    template_polytopes_m[a_type_id] = central_polytope_storage_m;
    needs_updating_m[a_type_id] = 0;
  }

  template <class GeometryType>
  PolytopeType &getHalfEdgePolytopeBase(const GeometryType &a_geometry) {
    const std::size_t type_id = this->getID(a_geometry);
    assert(type_id < needs_updating_m.size());
    assert(type_id < template_polytopes_m.size());
    if (needs_updating_m[type_id] == 1) {
      this->resetTemplatePolytope(a_geometry, type_id);
      return central_polytope_storage_m;
    }
    central_polytope_storage_m = template_polytopes_m[type_id];
    return central_polytope_storage_m;
  }

  ~HalfEdgeGeometryInitializer(void) = default;

private:
  PolytopeType central_polytope_storage_m;
  SmallVector<PolytopeType, 8> template_polytopes_m;
  SmallVector<UnsignedIndex_t, 8> needs_updating_m;
};

template <class VertexType>
HalfEdgeGeometryInitializer<HalfEdgePolyhedron<VertexType>> &
getHalfEdgePolyhedronStorage(void);

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            void>
updatePolytopeStorage(void);

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            void>
updatePolytopeStorage(void);

template <class EncompassingGeometryType>
enable_if_t<is_tri<EncompassingGeometryType>::value, void>
updatePolytopeStorage(void);

template <class EncompassingGeometryType>
enable_if_t<is_polygon<EncompassingGeometryType>::value &&
                !is_tri<EncompassingGeometryType>::value,
            void>
updatePolytopeStorage(void);

template <class VertexType, class GeometryType>
HalfEdgePolyhedron<VertexType> &getHalfEdgePolyhedron(void);

template <class VertexType>
HalfEdgePolygon<VertexType> &getHalfEdgePolygonStorage(void);

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry);

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry);

template <class EncompassingGeometryType>
enable_if_t<is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry);

template <class EncompassingGeometryType>
enable_if_t<is_polygon<EncompassingGeometryType>::value &&
                !is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry);

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolyhedron<typename PolytopeType::face_type,
                                        typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType *a_complete_polytope);

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polygon<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolygon<typename PolytopeType::face_type,
                                     typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType *a_complete_polytope);

//******************************************************************* //
//     Function template definitions placed below this
//******************************************************************* //
template <class ReturnType, class EncompassingType, class ReconstructionType>
ReturnType
cutThroughHalfEdgeStructures(const EncompassingType &a_polytope,
                             const ReconstructionType &a_reconstruction) {
  auto &complete_polytope = setHalfEdgeStructure(a_polytope);
  auto half_edge_polytope =
      generateSegmentedVersion<EncompassingType>(&complete_polytope);
  assert(half_edge_polytope.checkValidHalfEdgeStructure());
  const ReturnType moments = getVolumeMoments<ReturnType, HalfEdgeCutting>(
      &half_edge_polytope, &complete_polytope, a_reconstruction);
  updatePolytopeStorage<EncompassingType>();
  return moments;
}

template <class VertexType>
HalfEdgeGeometryInitializer<HalfEdgePolyhedron<VertexType>> &
getHalfEdgePolyhedronStorage(void) {
  static HalfEdgeGeometryInitializer<HalfEdgePolyhedron<VertexType>> storage;
  return storage;
}

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            void>
updatePolytopeStorage(void) {
  auto &storage = getHalfEdgePolyhedronStorage<
      typename EncompassingGeometryType::pt_type>();
  storage.resetCentralStorageToCurrentSize();
}

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            void>
updatePolytopeStorage(void) {}

template <class EncompassingGeometryType>
enable_if_t<is_tri<EncompassingGeometryType>::value, void>
updatePolytopeStorage(void) {}

template <class EncompassingGeometryType>
enable_if_t<is_polygon<EncompassingGeometryType>::value &&
                !is_tri<EncompassingGeometryType>::value,
            void>
updatePolytopeStorage(void) {}

template <class VertexType, class GeometryType>
HalfEdgePolyhedron<VertexType> &
getHalfEdgePolyhedron(const GeometryType &a_geometry) {
  auto &storage = getHalfEdgePolyhedronStorage<VertexType>();
  return storage.getHalfEdgePolytopeBase(a_geometry);
}

template <class VertexType>
HalfEdgePolygon<VertexType> &getHalfEdgePolygonStorage(void) {
  static HalfEdgePolygon<VertexType> storage;
  return storage;
}

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                !is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry) {
  auto &complete_polyhedron_buffer =
      getHalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>(
          a_geometry);
  complete_polyhedron_buffer.setVertexLocations(a_geometry);
  return complete_polyhedron_buffer;
}

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value &&
                is_general_polyhedron<EncompassingGeometryType>::value,
            HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry) {
  // Can't lazy evaluate GeneralPolyhedron's because they change, i.e.
  // since a GeneralPolyhedron doesn't have a fixed number of vertices or
  // connectivity like other polyhedrons in IRL, the connectivity cannot
  // simply be stored and copied over.
  static HalfEdgePolyhedron<typename EncompassingGeometryType::pt_type>
      complete_polyhedron_buffer;
  a_geometry.setHalfEdgeVersion(&complete_polyhedron_buffer);
  return complete_polyhedron_buffer;
}

template <class EncompassingGeometryType>
enable_if_t<is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry) {
  static bool already_set = false;
  static HalfEdgePolygon<typename EncompassingGeometryType::pt_type>
      half_edge_geometry_template;

  static HalfEdgePolygon<typename EncompassingGeometryType::pt_type>
      complete_polygon_buffer;

  if (already_set) {
    complete_polygon_buffer = half_edge_geometry_template;
    complete_polygon_buffer.setVertexLocations(a_geometry);
  } else {
    already_set = true;
    a_geometry.setHalfEdgeVersion(&complete_polygon_buffer);
    half_edge_geometry_template = complete_polygon_buffer;
  }
  complete_polygon_buffer.setPlaneOfExistence(a_geometry.getPlaneOfExistence());
  return complete_polygon_buffer;
}

template <class EncompassingGeometryType>
enable_if_t<is_polygon<EncompassingGeometryType>::value &&
                !is_tri<EncompassingGeometryType>::value,
            HalfEdgePolygon<typename EncompassingGeometryType::pt_type> &>
setHalfEdgeStructure(const EncompassingGeometryType &a_geometry) {
  // Can't lazy evaluate polygons because they change, i.e.
  // since a Polygon doesn't have a fixed number of vertices like
  // a polyhedron does in IRL, the connectivity cannot simply
  // be stored and copied over.
  HalfEdgePolygon<
      typename EncompassingGeometryType::pt_type> &complete_polygon_buffer =
      getHalfEdgePolygonStorage<typename EncompassingGeometryType::pt_type>();
  a_geometry.setHalfEdgeVersion(&complete_polygon_buffer);
  return complete_polygon_buffer;
}

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolyhedron<typename PolytopeType::face_type,
                                        typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType *a_complete_polytope) {
  return a_complete_polytope->generateSegmentedPolyhedron();
}

template <class EncompassingGeometryType, class PolytopeType>
enable_if_t<is_polygon<EncompassingGeometryType>::value,
            SegmentedHalfEdgePolygon<typename PolytopeType::face_type,
                                     typename PolytopeType::vertex_type>>
generateSegmentedVersion(PolytopeType *a_complete_polytope) {
  auto segmented_half_edge_geometry_template =
      a_complete_polytope->generateSegmentedPolygon();
  segmented_half_edge_geometry_template.setPlaneOfExistence(
      &a_complete_polytope->getPlaneOfExistence());
  return segmented_half_edge_geometry_template;
}

} // namespace IRL

#endif // SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_INITIALIZER_TPP_
