// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_CAPPED_DODECAHEDRON_LTLT_H_
#define IRL_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_CAPPED_DODECAHEDRON_LTLT_H_

#include "irl/geometry/general/stored_vertex_access.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/parameters/defined_types.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_correction_base.h"

namespace IRL{

template <class Derived, class VertexType>
class CappedDodecahedron_LTLTSpecialization: public CappedDodecahedronCorrectionBase<Derived, VertexType> {
public: 

HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const; 

template<class HalfEdgePolyhedronType>
void setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const;

static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

static constexpr std::array<UnsignedIndex_t, 4> getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

ProxyTet<Derived> getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const;
}; 

template <class VertexType>
class StoredCappedDodecahedron_LTLT: public StoredVertexAccess<StoredCappedDodecahedron_LTLT<VertexType>,VertexType,9>, public CappedDodecahedron_LTLTSpecialization<StoredCappedDodecahedron_LTLT<VertexType>, VertexType>{
friend StoredVertexAccess<StoredCappedDodecahedron_LTLT<VertexType>,VertexType,9>;

public: 

using StoredVertexAccess<StoredCappedDodecahedron_LTLT<VertexType>,VertexType,9>::StoredVertexAccess;

StoredCappedDodecahedron_LTLT(void) = default;
}; 


// Predefined types 
using CappedDodecahedron_LTLT= StoredCappedDodecahedron_LTLT<Pt>; 

} // namespace IRL 

#include "capped_dodecahedron_LTLT.tpp"
#endif //SRC_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_CAPPED_DODECAHEDRON_LTLT_H_
