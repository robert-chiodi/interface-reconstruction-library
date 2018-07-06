// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_DECOMPOSED_POLYTOPE_SEGMENTED_DECOMPOSED_POLYTOPE_H_
#define SRC_GEOMETRY_DECOMPOSED_POLYTOPE_SEGMENTED_DECOMPOSED_POLYTOPE_H_

#include <algorithm>
#include <numeric>
#include <utility>

#include "src/data_structures/stack_vector.h"
#include "src/geometry/decomposed_polytope/decomposed_polytope_vertex_storage.h"
#include "src/geometry/general/geometry_type_traits.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class ContainerType>
class IteratorThroughBracketOperator;

template <class ContainerType>
class ConstIteratorThroughBracketOperator;

template <class VertexStorage, class VertexType>
class VertexList {
 public:
  using iterator =
      IteratorThroughBracketOperator<VertexList<VertexStorage, VertexType>>;
  using const_iterator = ConstIteratorThroughBracketOperator<
      VertexList<VertexStorage, VertexType>>;
  using value_t = VertexType;

  explicit VertexList(VertexStorage& a_reference)
      : reference_list_m(a_reference) {}

  void assumeWholeOfReference(void) {
    referenced_vertices_m.resize(reference_list_m.size());
    std::iota(referenced_vertices_m.begin(),
              referenced_vertices_m.begin() + reference_list_m.size(), 0);
  }

  VertexList(const VertexList& a_other) noexcept
      : reference_list_m(a_other.reference_list_m) {
    referenced_vertices_m = a_other.referenced_vertices_m;
  }

  VertexList(VertexList&& a_other) noexcept
      : reference_list_m(a_other.reference_list_m) {
    referenced_vertices_m = std::move(a_other.referenced_vertices_m);
  }

  VertexList& operator=(const VertexList& a_other) noexcept {
    assert(&reference_list_m == &a_other.reference_list_m);
    if (this != &a_other) {
      referenced_vertices_m = a_other.referenced_vertices_m;
    }
    return (*this);
  }

  VertexList& operator=(VertexList&& a_other) noexcept {
    assert(&reference_list_m == &a_other.reference_list_m);
    if (this != &a_other) {
      referenced_vertices_m = std::move(a_other.referenced_vertices_m);
    }
    return (*this);
  }

  VertexType& operator[](const UnsignedIndex_t a_index) {
    return reference_list_m[referenced_vertices_m[a_index]];
  }
  const VertexType& operator[](const UnsignedIndex_t a_index) const {
    return reference_list_m[referenced_vertices_m[a_index]];
  }

  double getDistance(const UnsignedIndex_t a_index) const {
    return reference_list_m.getDistance(referenced_vertices_m[a_index]);
  }

  void setDistance(const UnsignedIndex_t a_index,
                   const double a_distance) const {
    reference_list_m.setDistance(referenced_vertices_m[a_index], a_distance);
  }

  UnsignedIndex_t& getVertexIndex(const UnsignedIndex_t a_index) {
    return referenced_vertices_m[a_index];
  }

  const UnsignedIndex_t& getVertexIndex(const UnsignedIndex_t a_index) const {
    return referenced_vertices_m[a_index];
  }

  void rePointVertex(const UnsignedIndex_t a_index,
                     const UnsignedIndex_t a_new_vertex_index) {
    assert(a_new_vertex_index < reference_list_m.size());
    referenced_vertices_m[a_index] = a_new_vertex_index;
  }

  void push_back(const UnsignedIndex_t a_new_vertex_index) {
    referenced_vertices_m.push_back(a_new_vertex_index);
  }

  UnsignedIndex_t size(void) const {
    return static_cast<UnsignedIndex_t>(referenced_vertices_m.size());
  }
  void resize(const UnsignedIndex_t a_size) {
    referenced_vertices_m.resize(a_size);
  }

  VertexStorage& getUnderlyingVertexStorage(void) { return reference_list_m; }
  const VertexStorage& getUnderlyingVertexStorage(void) const {
    return reference_list_m;
  }

  iterator begin(void) noexcept { return iterator((*this), 0); }
  const_iterator begin(void) const noexcept { return this->cbegin(); }
  const_iterator cbegin(void) const noexcept {
    return const_iterator((*this), 0);
  }
  iterator end(void) noexcept { return iterator((*this), this->size()); }
  const_iterator end(void) const noexcept { return this->cend(); }
  const_iterator cend(void) const noexcept {
    return const_iterator((*this), this->size());
  }

 private:
  VertexStorage& reference_list_m;
  SmallVector<UnsignedIndex_t, 40> referenced_vertices_m;
};

template <class VertexType>
class SegmentedDecomposedPolyhedron;

template <class VertexType>
class SegmentedDecomposedPolygon;

template <class VertexType, class ProxyType, typename Enable = void>
class MomentCalculationType;

template <class VertexType, class ProxyType>
class MomentCalculationType<VertexType, ProxyType,
                            enable_if_t<is_polyhedron<ProxyType>::value>>
    : public PolyhedronMomentsCalculation<
          SegmentedDecomposedPolyhedron<VertexType>, VertexType, ProxyType> {};

template <class VertexType, class ProxyType>
class MomentCalculationType<VertexType, ProxyType,
                            enable_if_t<is_polygon<ProxyType>::value>>
    : public PolygonMomentsCalculation<SegmentedDecomposedPolygon<VertexType>,
                                       VertexType, ProxyType> {};

template <class VertexStorageType, class VertexType, class ProxyType>
class SegmentedDecomposedPolytope
    : public MomentCalculationType<VertexType, ProxyType> {
 public:
  using pt_type = VertexType;

  explicit SegmentedDecomposedPolytope(VertexStorageType& a_vertex_storage)
      : vertex_list_m(a_vertex_storage) {}

  VertexList<VertexStorageType, VertexType>& getVertexList(void) {
    return vertex_list_m;
  }

  const VertexList<VertexStorageType, VertexType>& getVertexList(void) const {
    return vertex_list_m;
  }

  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const {
    return static_cast<UnsignedIndex_t>(simplex_decomposition_m.size());
  }

  ProxyType getNullProxyObject(void) const {
    std::array<UnsignedIndex_t, ProxyType::getNumberOfVerticesInObject()>
        indices;
    // Fill with largest UnsignedIndex value
    std::fill(indices.begin(), indices.end(), static_cast<UnsignedIndex_t>(-1));
    return ProxyType(this->getVertexList().getUnderlyingVertexStorage(),
                     indices);
  }

  void setNumberOfSimplicesInDecomposition(const UnsignedIndex_t a_new_size) {
    simplex_decomposition_m.resize(a_new_size, this->getNullProxyObject());
  }

  ProxyType& getSimplexFromDecomposition(const UnsignedIndex_t a_tet) {
    assert(a_tet < this->getNumberOfSimplicesInDecomposition());
    return simplex_decomposition_m[a_tet];
  }

  const ProxyType& getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet) const {
    assert(a_tet < this->getNumberOfSimplicesInDecomposition());
    return simplex_decomposition_m[a_tet];
  }

  void setNumberOfValidSimplices(const UnsignedIndex_t a_new_size) {
    valid_simplex_mask.resize(a_new_size);
  }

  void markAsValid(const UnsignedIndex_t a_index) {
    valid_simplex_mask[a_index] = 1;
  }

  void markAsNotValid(const UnsignedIndex_t a_index) {
    valid_simplex_mask[a_index] = 0;
  }

  bool isValid(const UnsignedIndex_t a_index) const {
    return valid_simplex_mask[a_index] == 1;
  }

  bool isNotValid(const UnsignedIndex_t a_index) const {
    return !this->isValid(a_index);
  }

  void clear(void) {
    this->setNumberOfSimplicesInDecomposition(0);
    this->setNumberOfValidSimplices(0);
    vertex_list_m.resize(0);
  }

  VertexType& operator[](const UnsignedIndex_t a_index) {
    return vertex_list_m[a_index];
  }
  const VertexType& operator[](const UnsignedIndex_t a_index) const {
    return vertex_list_m[a_index];
  }

  typename VertexList<VertexStorageType, VertexType>::iterator begin(
      void) noexcept {
    return vertex_list_m.begin();
  }
  typename VertexList<VertexStorageType, VertexType>::const_iterator begin(
      void) const noexcept {
    return this->cbegin();
  }
  typename VertexList<VertexStorageType, VertexType>::const_iterator cbegin(
      void) const noexcept {
    return vertex_list_m.cbegin();
  }
  typename VertexList<VertexStorageType, VertexType>::iterator end(
      void) noexcept {
    return vertex_list_m.end();
  }
  typename VertexList<VertexStorageType, VertexType>::const_iterator end(
      void) const noexcept {
    return this->cend();
  }
  typename VertexList<VertexStorageType, VertexType>::const_iterator cend(
      void) const noexcept {
    return vertex_list_m.cend();
  }

 protected:
  SmallVector<ProxyType, 40> simplex_decomposition_m;
  SmallVector<uint8_t, 40> valid_simplex_mask;
  VertexList<VertexStorageType, VertexType> vertex_list_m;
};

template <class VertexStorageType, class VertexType, class ProxyType>
inline std::ostream& operator<<(
    std::ostream& out,
    const SegmentedDecomposedPolytope<VertexStorageType, VertexType, ProxyType>&
        a_decomposed_polytope) {
  out << "Polytope object made of "
      << a_decomposed_polytope.getNumberOfSimplicesInDecomposition()
      << " simplices " << std::endl;
  out << "They are: " << std::endl;
  for (UnsignedIndex_t n = 0;
       n < a_decomposed_polytope.getNumberOfSimplicesInDecomposition(); ++n) {
    std::cout << "Simplex " << n << std::endl;
    std::cout << a_decomposed_polytope.getSimplexFromDecomposition(n) << '\n'
              << std::endl;
  }
  out << "=============================================\n" << std::endl;
  return out;
}

template <class VertexType>
class SegmentedDecomposedPolyhedron
    : public SegmentedDecomposedPolytope<
          DecomposedPolyhedronVertexStorage<VertexType>, VertexType,
          ProxyTet<DecomposedPolyhedronVertexStorage<VertexType>>> {
 public:
  using VertexStorageType = DecomposedPolyhedronVertexStorage<VertexType>;
  using ProxyType = ProxyTet<DecomposedPolyhedronVertexStorage<VertexType>>;
  using BaseClass = SegmentedDecomposedPolytope<
      DecomposedPolyhedronVertexStorage<VertexType>, VertexType,
      ProxyTet<DecomposedPolyhedronVertexStorage<VertexType>>>;

  template <class GeometryType>
  SegmentedDecomposedPolyhedron(const GeometryType& a_geometry,
                                VertexStorageType& a_vertex_storage)
      : BaseClass(a_vertex_storage) {
    this->setNumberOfSimplicesInDecomposition(
        a_geometry.getNumberOfSimplicesInDecomposition());
    for (UnsignedIndex_t n = 0;
         n < a_geometry.getNumberOfSimplicesInDecomposition(); ++n) {
      this->getSimplexFromDecomposition(n) = ProxyType(
          a_vertex_storage, a_geometry.getSimplexIndicesFromDecomposition(n));
    }
    this->vertex_list_m.assumeWholeOfReference();
    this->valid_simplex_mask.assign(
        a_geometry.getNumberOfSimplicesInDecomposition(), true);
  }

  explicit SegmentedDecomposedPolyhedron(VertexStorageType& a_vertex_storage)
      : BaseClass(a_vertex_storage) {}

  void push_back(const ProxyType& a_proxy) {
    this->simplex_decomposition_m.push_back(a_proxy);
    this->valid_simplex_mask.push_back(true);
  }

 private:
};

template <class VertexType>
class SegmentedDecomposedPolygon
    : public SegmentedDecomposedPolytope<
          DecomposedPolygonVertexStorage<VertexType>, VertexType,
          ProxyTri<DecomposedPolygonVertexStorage<VertexType>>> {
 public:
  using VertexStorageType = DecomposedPolygonVertexStorage<VertexType>;
  using ProxyType = ProxyTri<DecomposedPolygonVertexStorage<VertexType>>;
  using BaseClass = SegmentedDecomposedPolytope<
      DecomposedPolygonVertexStorage<VertexType>, VertexType,
      ProxyTri<DecomposedPolygonVertexStorage<VertexType>>>;

  template <class GeometryType>
  SegmentedDecomposedPolygon(const GeometryType& a_geometry,
                             VertexStorageType& a_vertex_storage)
      : BaseClass(a_vertex_storage),
        plane_of_existence_m(&a_vertex_storage.getPlaneOfExistence()) {
    this->setNumberOfSimplicesInDecomposition(
        a_geometry.getNumberOfSimplicesInDecomposition());
    for (UnsignedIndex_t n = 0;
         n < a_geometry.getNumberOfSimplicesInDecomposition(); ++n) {
      this->getSimplexFromDecomposition(n) = ProxyType(
          a_vertex_storage, a_geometry.getSimplexIndicesFromDecomposition(n),
          this->getPlaneOfExistence());
    }
    this->vertex_list_m.assumeWholeOfReference();
    this->valid_simplex_mask.assign(
        a_geometry.getNumberOfSimplicesInDecomposition(), true);
  }

  explicit SegmentedDecomposedPolygon(VertexStorageType& a_vertex_storage)
      : BaseClass(a_vertex_storage),
        plane_of_existence_m(&a_vertex_storage.getPlaneOfExistence()) {}

  ProxyType getNullProxyObject(void) const {
    std::array<UnsignedIndex_t, ProxyType::getNumberOfVerticesInObject()>
        indices;
    // Fill with largest UnsignedIndex value
    std::fill(indices.begin(), indices.end(), static_cast<UnsignedIndex_t>(-1));
    return ProxyType(this->getVertexList().getUnderlyingVertexStorage(),
                     indices, this->getPlaneOfExistence());
  }

  void setNumberOfSimplicesInDecomposition(const UnsignedIndex_t a_new_size) {
    this->simplex_decomposition_m.resize(a_new_size,
                                         this->getNullProxyObject());
  }

  void push_back(const ProxyType& a_proxy) {
    this->simplex_decomposition_m.push_back(a_proxy);
    this->valid_simplex_mask.push_back(true);
  }

  void setPlaneOfExistence(const Plane* a_plane) {
    plane_of_existence_m = a_plane;
  }
  const Plane& getPlaneOfExistence(void) const {
    assert(plane_of_existence_m != nullptr);
    return *plane_of_existence_m;
  }

  VolumeMomentsAndNormal calculateVolumeMomentsAndNormal(void) const {
    return calculateMoments((*this), VolumeMomentsAndNormal2D_Functor());
  }

 private:
  const Plane* plane_of_existence_m;
};

}  // namespace IRL

#endif  // SRC_GEOMETRY_DECOMPOSED_POLYTOPE_SEGMENTED_DECOMPOSED_POLYTOPE_H_
