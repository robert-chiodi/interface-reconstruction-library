// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_DATA_STRUCTURES_ITERATOR_THROUGH_BRACKET_OPERATOR_H_
#define IRL_DATA_STRUCTURES_ITERATOR_THROUGH_BRACKET_OPERATOR_H_

#include <iterator>

#include "irl/parameters/defined_types.h"

namespace IRL {

template <class ContainerType>
class IteratorThroughBracketOperator {
 public:
  using value_t = typename ContainerType::value_t;
  using iterator_category = std::random_access_iterator_tag;
  using value_type = typename ContainerType::value_t;
  using difference_type = std::ptrdiff_t;
  using pointer = value_t*;
  using reference = value_t&;

  IteratorThroughBracketOperator(ContainerType& a_container,
                                 const UnsignedIndex_t a_location);

  value_t& operator*(void);
  value_t& operator->(void);

  value_t& operator[](const std::ptrdiff_t a_index);

  IteratorThroughBracketOperator& operator++(void);
  IteratorThroughBracketOperator& operator--(void);

  IteratorThroughBracketOperator operator++(int a_dummy_for_postfix);
  IteratorThroughBracketOperator operator--(int a_dummy_for_postfix);

  IteratorThroughBracketOperator& operator+=(const std::ptrdiff_t a_shift);
  IteratorThroughBracketOperator& operator-=(const std::ptrdiff_t a_shift);

  std::ptrdiff_t operator-(const IteratorThroughBracketOperator& a_rhs);

  bool operator==(const IteratorThroughBracketOperator& a_rhs) const;
  bool operator!=(const IteratorThroughBracketOperator& a_rhs) const;
  bool operator<(const IteratorThroughBracketOperator& a_rhs) const;
  bool operator>(const IteratorThroughBracketOperator& a_rhs) const;
  bool operator>=(const IteratorThroughBracketOperator& a_rhs) const;
  bool operator<=(const IteratorThroughBracketOperator& a_rhs) const;

 private:
  ContainerType& container_m;
  UnsignedIndex_t location_m;
};

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator+(
    const IteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment);

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator+(
    const std::ptrdiff_t a_increment,
    const IteratorThroughBracketOperator<ContainerType>& a_iterator);

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator-(
    const IteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment);

template <class ContainerType>
IteratorThroughBracketOperator<ContainerType> operator-(
    const std::ptrdiff_t a_increment,
    const IteratorThroughBracketOperator<ContainerType>& a_iterator);

template <class ContainerType>
class ConstIteratorThroughBracketOperator {
 public:
  using value_t = typename ContainerType::value_t;
  using iterator_category = std::random_access_iterator_tag;
  using value_type = typename ContainerType::value_t;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_t*;
  using reference = const value_t&;

  ConstIteratorThroughBracketOperator(const ContainerType& a_container,
                                      const UnsignedIndex_t a_location);
  const value_t& operator*(void);
  const value_t& operator->(void);

  const value_t& operator[](const std::ptrdiff_t a_index);

  ConstIteratorThroughBracketOperator& operator++(void);
  ConstIteratorThroughBracketOperator& operator--(void);

  ConstIteratorThroughBracketOperator operator++(int a_dummy_for_postfix);
  ConstIteratorThroughBracketOperator operator--(int a_dummy_for_postfix);

  ConstIteratorThroughBracketOperator& operator+=(const std::ptrdiff_t a_shift);
  ConstIteratorThroughBracketOperator& operator-=(const std::ptrdiff_t a_shift);

  std::ptrdiff_t operator-(const ConstIteratorThroughBracketOperator& a_rhs);

  bool operator==(const ConstIteratorThroughBracketOperator& a_rhs) const;
  bool operator!=(const ConstIteratorThroughBracketOperator& a_rhs) const;
  bool operator<(const ConstIteratorThroughBracketOperator& a_rhs) const;
  bool operator>(const ConstIteratorThroughBracketOperator& a_rhs) const;
  bool operator>=(const ConstIteratorThroughBracketOperator& a_rhs) const;
  bool operator<=(const ConstIteratorThroughBracketOperator& a_rhs) const;

 private:
  const ContainerType& container_m;
  UnsignedIndex_t location_m;
};

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator+(
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment);

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator+(
    const std::ptrdiff_t a_increment,
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator);

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator-(
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator,
    const std::ptrdiff_t a_increment);

template <class ContainerType>
ConstIteratorThroughBracketOperator<ContainerType> operator-(
    const std::ptrdiff_t a_increment,
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator);

template <class ContainerType>
std::ptrdiff_t operator-(
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator_0,
    const ConstIteratorThroughBracketOperator<ContainerType>& a_iterator_1);

}  // namespace IRL

#include "irl/data_structures/iterator_through_bracket_operator.tpp"

#endif  // IRL_DATA_STRUCTURES_ITERATOR_THROUGH_BRACKET_OPERATOR_H_
