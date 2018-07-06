// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_RAW_POINTER_WRAPPER_H_
#define SRC_DATA_STRUCTURES_RAW_POINTER_WRAPPER_H_

#include "src/parameters/defined_types.h"

namespace IRL {

template<class ObjectType>
class RawPointerWrapper;

template<class ObjectType>
class ConstRawPointerWrapper;

/// \brief Helper factory function to auto determine correct type.
template<class ObjectType>
RawPointerWrapper<ObjectType> wrapRawPointer(ObjectType* a_object_ptr, const UnsignedIndex_t a_size);

/// \brief Helper factory function to auto determine correct type.
template<class ObjectType>
ConstRawPointerWrapper<ObjectType> wrapRawPointer(const ObjectType* a_object_ptr, const UnsignedIndex_t a_size);

/// \brief This class wraps a pointer with its
/// size and provides access through the []
/// operator and the random access iterator.
/// Mostly used when receiving information
/// from C/Fortran.
template <class ObjectType>
class RawPointerWrapper {
 public:
  using iterator = ObjectType*;
  using const_iterator = const ObjectType*;
  /// \brief Default constructor.
  RawPointerWrapper(void) = delete;

  RawPointerWrapper(ObjectType* a_object_ptr, const UnsignedIndex_t a_size);

  ObjectType& operator[](const UnsignedIndex_t a_index);
  const ObjectType& operator[](const UnsignedIndex_t a_index) const;
  UnsignedIndex_t size(void) const;

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~RawPointerWrapper(void) = default;

 private:
  ObjectType* ptr_m;
  UnsignedIndex_t size_m;
};

template <class ObjectType>
class ConstRawPointerWrapper {
 public:
  using const_iterator = const ObjectType*;
  /// \brief Default constructor.
  ConstRawPointerWrapper(void) = delete;

  ConstRawPointerWrapper(const ObjectType* a_object_ptr, const UnsignedIndex_t a_size);

  const ObjectType& operator[](const UnsignedIndex_t a_index) const;
  UnsignedIndex_t size(void) const;

  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~ConstRawPointerWrapper(void) = default;

 private:
  const ObjectType* ptr_m;
  UnsignedIndex_t size_m;
};

}  // namespace IRL

#include "src/data_structures/raw_pointer_wrapper.tpp"

#endif  // SRC_DATA_STRUCTURES_RAW_POINTER_WRAPPER_H_
