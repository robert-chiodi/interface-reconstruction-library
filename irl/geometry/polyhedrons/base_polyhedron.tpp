// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_BASE_POLYHEDRON_TPP_
#define IRL_GEOMETRY_POLYHEDRONS_BASE_POLYHEDRON_TPP_

namespace IRL {
template <class Derived, class VertexType, class SimplexType>
Derived& BasePolyhedron<Derived, VertexType, SimplexType>::getDerived(void) {
  return static_cast<Derived&>(*this);
}
template <class Derived, class VertexType, class SimplexType>
const Derived& BasePolyhedron<Derived, VertexType, SimplexType>::getDerived(
    void) const {
  return static_cast<const Derived&>(*this);
}
template <class Derived, class VertexType, class SimplexType>
template <class OtherPolytope>
Derived BasePolyhedron<Derived, VertexType, SimplexType>::fromOtherPolytope(
    const OtherPolytope& a_other_polytope) {
  Derived object_to_return;
  assert(object_to_return.getNumberOfVertices() ==
         a_other_polytope.getNumberOfVertices());
  for (UnsignedIndex_t n = 0; n < object_to_return.getNumberOfVertices(); ++n) {
    object_to_return[n] = a_other_polytope[n];
  }
  return object_to_return;
}
template <class Derived, class VertexType, class SimplexType>
VertexType& BasePolyhedron<Derived, VertexType, SimplexType>::operator[](
    const UnsignedIndex_t a_index) {
  return this->getDerived().access(a_index);
}

template <class Derived, class VertexType, class SimplexType>
const VertexType& BasePolyhedron<Derived, VertexType, SimplexType>::
operator[](const UnsignedIndex_t a_index) const {
  return this->getDerived().access(a_index);
}

template <class Derived, class VertexType, class SimplexType>
UnsignedIndex_t BasePolyhedron<Derived, VertexType, SimplexType>::
    getNumberOfSimplicesInDecomposition(void) const {
  return this->getDerived().getNumberOfSimplicesInDecomposition();
}

template <class Derived, class VertexType, class SimplexType>
constexpr std::array<UnsignedIndex_t, 4>
BasePolyhedron<Derived, VertexType, SimplexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  return Derived::getSimplexIndicesFromDecomposition(a_tet);
}

template <class Derived, class VertexType, class SimplexType>
SimplexType BasePolyhedron<Derived, VertexType, SimplexType>::
    getSimplexFromDecomposition(
        const UnsignedIndex_t a_tet_number_to_get) const {
  return this->getDerived().getSimplexFromDecomposition(a_tet_number_to_get);
}

template <class Derived, class VertexType, class SimplexType>
UnsignedIndex_t
BasePolyhedron<Derived, VertexType, SimplexType>::getNumberOfVertices(
    void) const {
  return this->getDerived().getNumberOfVerticesInObject();
}

template <class Derived, class VertexType, class SimplexType>
IRL::Pt BasePolyhedron<Derived, VertexType, SimplexType>::getLowerLimits(
    void) const {
  Pt pt_to_return(DBL_MAX, DBL_MAX, DBL_MAX);
  for (const auto& pt : (*this)) {
    pt_to_return[0] = std::min(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::min(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::min(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class Derived, class VertexType, class SimplexType>
IRL::Pt BasePolyhedron<Derived, VertexType, SimplexType>::getUpperLimits(
    void) const {
  Pt pt_to_return(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (const auto& pt : (*this)) {
    pt_to_return[0] = std::max(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::max(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::max(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class Derived, class VertexType, class SimplexType>
void BasePolyhedron<Derived, VertexType, SimplexType>::shift(
    const double a_x_shift, const double a_y_shift, const double a_z_shift) {
  Pt shift_by_pt(a_x_shift, a_y_shift, a_z_shift);
  for (auto& vertex : (*this)) {
    vertex += shift_by_pt;
  }
}

template <class Derived, class VertexType, class SimplexType>
typename BasePolyhedron<Derived, VertexType, SimplexType>::iterator
BasePolyhedron<Derived, VertexType, SimplexType>::begin(void) noexcept {
  return iterator(this->getDerived(), 0);
}
template <class Derived, class VertexType, class SimplexType>
typename BasePolyhedron<Derived, VertexType, SimplexType>::const_iterator
BasePolyhedron<Derived, VertexType, SimplexType>::begin(void) const
    noexcept {
  return this->cbegin();
}
template <class Derived, class VertexType, class SimplexType>
typename BasePolyhedron<Derived, VertexType, SimplexType>::const_iterator
BasePolyhedron<Derived, VertexType, SimplexType>::cbegin(void) const
    noexcept {
  return const_iterator(this->getDerived(), 0);
}
template <class Derived, class VertexType, class SimplexType>
typename BasePolyhedron<Derived, VertexType, SimplexType>::iterator
BasePolyhedron<Derived, VertexType, SimplexType>::end(void) noexcept {
  return iterator(this->getDerived(), this->getNumberOfVertices());
}
template <class Derived, class VertexType, class SimplexType>
typename BasePolyhedron<Derived, VertexType, SimplexType>::const_iterator
BasePolyhedron<Derived, VertexType, SimplexType>::end(void) const noexcept {
  return this->cend();
}
template <class Derived, class VertexType, class SimplexType>
typename BasePolyhedron<Derived, VertexType, SimplexType>::const_iterator
BasePolyhedron<Derived, VertexType, SimplexType>::cend(void) const noexcept {
  return const_iterator(this->getDerived(), this->getNumberOfVertices());
}

template <class Derived, class VertexType, class SimplexType>
std::ostream& operator<<(
    std::ostream& out,
    const BasePolyhedron<Derived, VertexType, SimplexType>& a_polyhedron) {
  out << "Polyhedron has " << a_polyhedron.getNumberOfVertices()
      << "vertices \n";
  for (UnsignedIndex_t vertex = 0; vertex < a_polyhedron.getNumberOfVertices();
       ++vertex) {
    out << "Vertex " << vertex << " : " << a_polyhedron[vertex] << '\n';
  }
  return out;
}

}  // namespace IRL

#endif // IRL_GEOMETRY_POLYHEDRONS_BASE_POLYHEDRON_TPP_
