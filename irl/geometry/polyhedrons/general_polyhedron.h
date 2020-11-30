// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_GENERAL_POLYHEDRON_H_
#define IRL_GEOMETRY_POLYHEDRONS_GENERAL_POLYHEDRON_H_

#include "irl/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "irl/geometry/polyhedrons/base_polyhedron.h"
#include "irl/geometry/polyhedrons/polyhedron_connectivity.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class Derived, class VertexType>
class GeneralPolyhedronSpecialization
    : public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
  friend Derived;

 public:
  GeneralPolyhedronSpecialization(void);

  GeneralPolyhedronSpecialization(const PolyhedronConnectivity* a_connectivity);

  HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const;

  template <class HalfEdgePolyhedronType>
  void setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const;

  std::array<UnsignedIndex_t, 4> getSimplexIndicesFromDecomposition(
      const UnsignedIndex_t a_tet) const;

  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const;

  ProxyTet<Derived> getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet) const;

  void resetConnectivity(const PolyhedronConnectivity* a_connectivity);

 private:
  const PolyhedronConnectivity* connectivity_m;
};

template <class VertexType, UnsignedIndex_t kStaticAllocSize = 10>
class StoredGeneralPolyhedron
    : public GeneralPolyhedronSpecialization<
          StoredGeneralPolyhedron<VertexType, kStaticAllocSize>, VertexType> {
  using Base = GeneralPolyhedronSpecialization<
      StoredGeneralPolyhedron<VertexType, kStaticAllocSize>, VertexType>;

 public:
  StoredGeneralPolyhedron(void) = default;

  template <class VertexListType>
  StoredGeneralPolyhedron(const VertexListType& a_vertex_list,
                          const PolyhedronConnectivity* a_connectivity);

  static StoredGeneralPolyhedron fromRawPtPointer(
      const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts,
      const PolyhedronConnectivity* a_connectivity);

  static StoredGeneralPolyhedron fromRawDoublePointer(
      const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs,
      const PolyhedronConnectivity* a_connectivity);

  VertexType& access(const UnsignedIndex_t a_index);
  const VertexType& access(const UnsignedIndex_t a_index) const;

  UnsignedIndex_t getNumberOfVerticesInObject(void) const;
  void setNumberOfVertices(const UnsignedIndex_t a_number);

 private:
  StoredGeneralPolyhedron(const UnsignedIndex_t a_number_of_pts,
                          const VertexType* a_array_of_pts,
                          const PolyhedronConnectivity* a_connectivity);

  StoredGeneralPolyhedron(const UnsignedIndex_t a_number_of_pts,
                          const double* a_array_of_locs,
                          const PolyhedronConnectivity* a_connectivity);

  ExpandablePtList<VertexType, kStaticAllocSize> vertex_list_m;
};

// Predefined types
using GeneralPolyhedron = StoredGeneralPolyhedron<Pt>;

}  // namespace IRL

#include "irl/geometry/polyhedrons/general_polyhedron.tpp"
#endif // IRL_GEOMETRY_POLYHEDRONS_GENERAL_POLYHEDRON_H_
