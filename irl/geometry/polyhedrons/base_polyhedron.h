// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_BASE_POLYHEDRON_H_
#define IRL_GEOMETRY_POLYHEDRONS_BASE_POLYHEDRON_H_

#include <algorithm>
#include <cstring>
#include <initializer_list>

#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/polyhedron_moments_calculation.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class ContainerType>
class IteratorThroughBracketOperator;

template <class ContainerType>
class ConstIteratorThroughBracketOperator;

template <class GeometryType>
class ProxyTet;

template <class Derived, class VertexType, class SimplexType>
class BasePolyhedron
    : public PolyhedronMomentsCalculation<Derived, VertexType, SimplexType> {
  Derived& getDerived(void);
  const Derived& getDerived(void) const;

  using iterator = IteratorThroughBracketOperator<Derived>;
  using const_iterator = ConstIteratorThroughBracketOperator<Derived>;

 public:
  using pt_type = VertexType;
  using value_t = pt_type;

  template <class OtherPolytope>
  static Derived fromOtherPolytope(const OtherPolytope& a_other_polytope);

  VertexType& operator[](const UnsignedIndex_t a_index);

  const VertexType& operator[](const UnsignedIndex_t a_index) const;

  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const;

  static constexpr std::array<UnsignedIndex_t, 4>
  getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

  SimplexType getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet_number_to_get) const;

  UnsignedIndex_t getNumberOfVertices(void) const;

  /// \brief Return a point for the lower limits of the polyhedron in 3D space.
  IRL::Pt getLowerLimits(void) const;

  /// \brief Return a point for the upper limits of the polyhedron in 3D space.
  IRL::Pt getUpperLimits(void) const;

  /// \brief Shift entire rectangular cuboid in x, y, and z.
  void shift(const double a_x_shift, const double a_y_shift,
             const double a_z_shift);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;
};

template <class Derived, class VertexType, class SimplexType>
std::ostream& operator<<(
    std::ostream& out,
    const BasePolyhedron<Derived, VertexType, SimplexType>& a_polyhedron);

}  // namespace IRL

#include "irl/geometry/polyhedrons/base_polyhedron.tpp"

#endif // IRL_GEOMETRY_POLYHEDRONS_BASE_POLYHEDRON_H_
