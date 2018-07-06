// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_RAW_POINTER_WRAPPER_TPP_
#define SRC_DATA_STRUCTURES_RAW_POINTER_WRAPPER_TPP_

#include <cassert>

namespace IRL {

template<class ObjectType>
RawPointerWrapper<ObjectType> wrapRawPointer(ObjectType* a_object_ptr, const UnsignedIndex_t a_size){
	return {a_object_ptr, a_size};
}

template<class ObjectType>
ConstRawPointerWrapper<ObjectType> wrapRawPointer(const ObjectType* a_object_ptr, const UnsignedIndex_t a_size){
	return {a_object_ptr, a_size};
}

template <class ObjectType>
RawPointerWrapper<ObjectType>::RawPointerWrapper(ObjectType* a_object_ptr,
												 const UnsignedIndex_t a_size) : ptr_m(a_object_ptr), size_m(a_size){
	assert(ptr_m != nullptr);
}

template <class ObjectType>
ObjectType& RawPointerWrapper<ObjectType>::operator[](const UnsignedIndex_t a_index){
	assert(a_index < size_m);
	return ptr_m[a_index];
}

template <class ObjectType>
const ObjectType& RawPointerWrapper<ObjectType>::operator[](const UnsignedIndex_t a_index) const{
	assert(a_index < size_m);
	return ptr_m[a_index];
}
template <class ObjectType>
UnsignedIndex_t RawPointerWrapper<ObjectType>::size(void) const{
  return size_m;
}

template <class ObjectType>
typename RawPointerWrapper<ObjectType>::iterator RawPointerWrapper<ObjectType>::begin(void) noexcept {
 return ptr_m;
}
template <class ObjectType>
typename RawPointerWrapper<ObjectType>::const_iterator RawPointerWrapper<ObjectType>::begin(void) const noexcept {
 return this->cbegin();
}
template <class ObjectType>
typename RawPointerWrapper<ObjectType>::const_iterator RawPointerWrapper<ObjectType>::cbegin(void) const noexcept {
 return ptr_m;
}

template <class ObjectType>
typename RawPointerWrapper<ObjectType>::iterator RawPointerWrapper<ObjectType>::end(void) noexcept {
 return ptr_m+size_m;
}
template <class ObjectType>
typename RawPointerWrapper<ObjectType>::const_iterator RawPointerWrapper<ObjectType>::end(void) const noexcept {
 return this->cend();
}
template <class ObjectType>
typename RawPointerWrapper<ObjectType>::const_iterator RawPointerWrapper<ObjectType>::cend(void) const noexcept {
 return ptr_m+size_m;
}

template <class ObjectType>
ConstRawPointerWrapper<ObjectType>::ConstRawPointerWrapper(const ObjectType* a_object_ptr,
												 const UnsignedIndex_t a_size) : ptr_m(a_object_ptr), size_m(a_size){
	assert(ptr_m != nullptr);
}


template <class ObjectType>
const ObjectType& ConstRawPointerWrapper<ObjectType>::operator[](const UnsignedIndex_t a_index) const{
	assert(a_index < size_m);
	return ptr_m[a_index];
}

template <class ObjectType>
UnsignedIndex_t ConstRawPointerWrapper<ObjectType>::size(void) const{
  return size_m;
}

template <class ObjectType>
typename ConstRawPointerWrapper<ObjectType>::const_iterator ConstRawPointerWrapper<ObjectType>::begin(void) const noexcept {
 return this->cbegin();
}
template <class ObjectType>
typename ConstRawPointerWrapper<ObjectType>::const_iterator ConstRawPointerWrapper<ObjectType>::cbegin(void) const noexcept {
 return ptr_m;
}

template <class ObjectType>
typename ConstRawPointerWrapper<ObjectType>::const_iterator ConstRawPointerWrapper<ObjectType>::end(void) const noexcept {
 return this->cend();
}
template <class ObjectType>
typename ConstRawPointerWrapper<ObjectType>::const_iterator ConstRawPointerWrapper<ObjectType>::cend(void) const noexcept {
 return ptr_m+size_m;
}


}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_RAW_POINTER_WRAPPER_TPP_
