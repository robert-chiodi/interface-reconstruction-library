// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_24_H_
#define SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_24_H_


#include "src/geometry/general/stored_vertex_access.h"
#include "src/parameters/defined_types.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_hexahedron.h"

namespace IRL {

template <class Derived, class VertexType>  
using Polyhedron24Specialization = SymmetricHexahedronSpecialization<Derived, VertexType>;

template <class VertexType>
class StoredPolyhedron24
    : public StoredVertexAccess<StoredPolyhedron24<VertexType>, VertexType, 14>,
      public Polyhedron24Specialization<StoredPolyhedron24<VertexType>,
                                        VertexType> {
  friend StoredVertexAccess<StoredPolyhedron24<VertexType>, VertexType, 14>;

 public:
  using StoredVertexAccess<StoredPolyhedron24<VertexType>, VertexType,
                           14>::StoredVertexAccess;

  StoredPolyhedron24(void) = default;
};

// Predefined types
using Polyhedron24 = StoredPolyhedron24<Pt>;

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_24_H_
