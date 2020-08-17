// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/planar_reconstruction/c_planar_separator_path_group.h"

#include <cassert>

#include "src/parameters/defined_types.h"

extern "C" {

void c_PlanarSepPathGroup_new(c_PlanarSepPathGroup* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PlanarSeparatorPathGroup;
}


void c_PlanarSepPathGroup_delete(c_PlanarSepPathGroup* a_self) {
  a_self->obj_ptr = nullptr;
}

void c_PlanarSepPathGroup_addPlanarSeparatorPath(c_PlanarSepPathGroup* a_self,
												 const c_PlanarSepPath* a_planar_separator_path){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_planar_separator_path != nullptr);
  assert(a_planar_separator_path->obj_ptr != nullptr);
  a_self->obj_ptr->addPlanarSeparatorPath(*(a_planar_separator_path->obj_ptr));
}

void c_PlanarSepPathGroup_addPlanarSeparatorPath_Id(c_PlanarSepPathGroup* a_self,
												 const c_PlanarSepPath* a_planar_separator_path,
												 const int* a_id){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_planar_separator_path != nullptr);
  assert(a_planar_separator_path->obj_ptr != nullptr);
  assert(*a_id >= 0);
  a_self->obj_ptr->addPlanarSeparatorPath(*(a_planar_separator_path->obj_ptr),
		  	  	  	  	  	  	  	      static_cast<IRL::UnsignedIndex_t>(*a_id));
}

void c_PlanarSepPathGroup_setPriorityOrder(c_PlanarSepPathGroup* a_self,
										   const int* a_size,
										   const int* a_priority_order){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_size >= 0);
  a_self->obj_ptr->setPriorityOrder(static_cast<IRL::UnsignedIndex_t>(*a_size), a_priority_order);
}

int c_PlanarSepPathGroup_getPriorityOrderSize(c_PlanarSepPathGroup* a_self){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getPriorityOrderSize());
}

int c_PlanarSepPathGroup_getPriorityOrderTag(c_PlanarSepPathGroup* a_self, const int* a_index){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getPriorityOrderTag(static_cast<IRL::UnsignedIndex_t>(*a_index)));
}

} // end extern "C"
