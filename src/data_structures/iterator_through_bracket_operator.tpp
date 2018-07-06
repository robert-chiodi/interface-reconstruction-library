// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_ITERATOR_THROUGH_BRACKET_OPERATOR_TPP_
#define SRC_DATA_STRUCTURES_ITERATOR_THROUGH_BRACKET_OPERATOR_TPP_

namespace IRL {

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType>::IteratorThroughBracketOperator(
    ContainerType& a_container, const UnsignedIndex_t a_location)
    : container_m{a_container}, location_m{a_location} {}

template <class ContainerType>
typename IteratorThroughBracketOperator<ContainerType>::value_t&
    IteratorThroughBracketOperator<ContainerType>::operator*(void) {
  return container_m[location_m];
}
template <class ContainerType>
typename IteratorThroughBracketOperator<ContainerType>::value_t&
    IteratorThroughBracketOperator<ContainerType>::operator->(void) {
  return container_m[location_m];
}

template <class ContainerType>
typename IteratorThroughBracketOperator<ContainerType>::value_t&
    IteratorThroughBracketOperator<ContainerType>::operator[](
        const std::ptrdiff_t a_index) {
  return container_m[location_m + a_index];
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType>&
IteratorThroughBracketOperator<ContainerType>::operator++(void) {
  ++location_m;
  return (*this);
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType>&
IteratorThroughBracketOperator<ContainerType>::operator--(void) {
  --location_m;
  return (*this);
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType>
IteratorThroughBracketOperator<ContainerType>::operator++(
    int a_dummy_for_postfix) {
  IteratorThroughBracketOperator copy_of_iterator(*this);
  ++location_m;
  return copy_of_iterator;
}
template <class ContainerType>
IteratorThroughBracketOperator<ContainerType>
IteratorThroughBracketOperator<ContainerType>::operator--(
    int a_dummy_for_postfix) {
  IteratorThroughBracketOperator copy_of_iterator(*this);
  --location_m;
  return copy_of_iterator;
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType>&
IteratorThroughBracketOperator<ContainerType>::operator+=(
    const std::ptrdiff_t a_shift) {
  location_m = static_cast<UnsignedIndex_t>(
      static_cast<std::ptrdiff_t>(location_m) + a_shift);
  return (*this);
}
template <class ContainerType>
IteratorThroughBracketOperator<ContainerType>&
IteratorThroughBracketOperator<ContainerType>::operator-=(
    const std::ptrdiff_t a_shift) {
  location_m = static_cast<UnsignedIndex_t>(
      static_cast<std::ptrdiff_t>(location_m) - a_shift);
  return (*this);
}

template <class ContainerType>
std::ptrdiff_t IteratorThroughBracketOperator<ContainerType>::operator-(
    const IteratorThroughBracketOperator& a_rhs) {
  return location_m - a_rhs.location_m;
}

template <class ContainerType>
bool IteratorThroughBracketOperator<ContainerType>::operator==(
    const IteratorThroughBracketOperator& a_rhs) const {
  return &container_m == &a_rhs.container_m && location_m == a_rhs.location_m;
}

template <class ContainerType>
bool IteratorThroughBracketOperator<ContainerType>::operator!=(
    const IteratorThroughBracketOperator& a_rhs) const {
  return !((*this) == a_rhs);
}
template <class ContainerType>
bool IteratorThroughBracketOperator<ContainerType>::operator<(
    const IteratorThroughBracketOperator& a_rhs) const {
  return location_m < a_rhs.location_m;
}
template <class ContainerType>
bool IteratorThroughBracketOperator<ContainerType>::operator>(
    const IteratorThroughBracketOperator& a_rhs) const {
  return a_rhs < (*this);
}
template <class ContainerType>
bool IteratorThroughBracketOperator<ContainerType>::operator>=(
    const IteratorThroughBracketOperator& a_rhs) const {
  return !((*this) < a_rhs);
}
template <class ContainerType>
bool IteratorThroughBracketOperator<ContainerType>::operator<=(
    const IteratorThroughBracketOperator& a_rhs) const {
  return !(a_rhs < (*this));
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator+(
    const IteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment) {
  IteratorThroughBracketOperator<ContainerType> copy = a_iterator;
  copy += a_increment;
  return copy;
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator+(
    const std::ptrdiff_t a_increment,
    const IteratorThroughBracketOperator<ContainerType>& a_iterator) {
  return a_iterator + a_increment;
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator-(
    const IteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment) {
  return a_iterator + (-a_increment);
}

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator-(
    const std::ptrdiff_t a_increment,
    const IteratorThroughBracketOperator<ContainerType>& a_iterator) {
  return a_iterator - a_increment;
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>::
    ConstIteratorThroughBracketOperator(const ContainerType& a_container,
                                        const UnsignedIndex_t a_location)
    : container_m{a_container}, location_m{a_location} {}

template <class ContainerType>
const typename ConstIteratorThroughBracketOperator<ContainerType>::value_t&
    ConstIteratorThroughBracketOperator<ContainerType>::operator*(void) {
  return container_m[location_m];
}
template <class ContainerType>
const typename ConstIteratorThroughBracketOperator<ContainerType>::value_t&
    ConstIteratorThroughBracketOperator<ContainerType>::operator->(void) {
  return container_m[location_m];
}

template <class ContainerType>
const typename ConstIteratorThroughBracketOperator<ContainerType>::value_t&
    ConstIteratorThroughBracketOperator<ContainerType>::operator[](
        const std::ptrdiff_t a_index) {
  return container_m[location_m + a_index];
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>&
ConstIteratorThroughBracketOperator<ContainerType>::operator++(void) {
  ++location_m;
  return (*this);
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>&
ConstIteratorThroughBracketOperator<ContainerType>::operator--(void) {
  --location_m;
  return (*this);
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>::operator++(
    int a_dummy_for_postfix) {
  ConstIteratorThroughBracketOperator copy_of_iterator(*this);
  ++location_m;
  return copy_of_iterator;
}
template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>::operator--(
    int a_dummy_for_postfix) {
  ConstIteratorThroughBracketOperator copy_of_iterator(*this);
  --location_m;
  return copy_of_iterator;
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>&
ConstIteratorThroughBracketOperator<ContainerType>::operator+=(
    const std::ptrdiff_t a_shift) {
  location_m = static_cast<UnsignedIndex_t>(
      static_cast<std::ptrdiff_t>(location_m) + a_shift);
  return (*this);
}
template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType>&
ConstIteratorThroughBracketOperator<ContainerType>::operator-=(
    const std::ptrdiff_t a_shift) {
  location_m = static_cast<UnsignedIndex_t>(
      static_cast<std::ptrdiff_t>(location_m) - a_shift);
  return (*this);
}

template <class ContainerType>
std::ptrdiff_t ConstIteratorThroughBracketOperator<ContainerType>::operator-(
    const ConstIteratorThroughBracketOperator& a_rhs) {
  return location_m - a_rhs.location_m;
}

template <class ContainerType>
bool ConstIteratorThroughBracketOperator<ContainerType>::operator==(
    const ConstIteratorThroughBracketOperator& a_rhs) const {
  return &container_m == &a_rhs.container_m && location_m == a_rhs.location_m;
}

template <class ContainerType>
bool ConstIteratorThroughBracketOperator<ContainerType>::operator!=(
    const ConstIteratorThroughBracketOperator& a_rhs) const {
  return !((*this) == a_rhs);
}
template <class ContainerType>
bool ConstIteratorThroughBracketOperator<ContainerType>::operator<(
    const ConstIteratorThroughBracketOperator& a_rhs) const {
  return location_m < a_rhs.location_m;
}
template <class ContainerType>
bool ConstIteratorThroughBracketOperator<ContainerType>::operator>(
    const ConstIteratorThroughBracketOperator& a_rhs) const {
  return a_rhs < (*this);
}
template <class ContainerType>
bool ConstIteratorThroughBracketOperator<ContainerType>::operator>=(
    const ConstIteratorThroughBracketOperator& a_rhs) const {
  return !((*this) < a_rhs);
}
template <class ContainerType>
bool ConstIteratorThroughBracketOperator<ContainerType>::operator<=(
    const ConstIteratorThroughBracketOperator& a_rhs) const {
  return !(a_rhs < (*this));
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator+(
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment) {
  ConstIteratorThroughBracketOperator<ContainerType> copy = a_iterator;
  copy += a_increment;
  return copy;
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator+(
    const std::ptrdiff_t a_increment,
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator) {
  return a_iterator + a_increment;
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator-(
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment) {
  return a_iterator + (-a_increment);
}

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator-(
    const std::ptrdiff_t a_increment,
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator) {
  return a_iterator - a_increment;
}

template <class ContainerType>
std::ptrdiff_t operator-(
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator_0,
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator_1) {
  return a_iterator_0.location_m - a_iterator_1.location_m;
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_ITERATOR_THROUGH_BRACKET_OPERATOR_TPP_
