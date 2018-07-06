// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_GROUP_TPP_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_GROUP_TPP_

#include <cassert>
#include <utility>

namespace IRL {

inline void PlanarSeparatorPathGroup::addPlanarSeparatorPath(
    const PlanarSeparatorPath& a_planar_separator_path) {
  assert(a_planar_separator_path.isIdSet());
  assert(this->isTagNew(a_planar_separator_path.getId()));
  id_mapping_m.insert(
      {a_planar_separator_path.getId(), a_planar_separator_path});
}
inline void PlanarSeparatorPathGroup::addPlanarSeparatorPath(
    PlanarSeparatorPath&& a_planar_separator_path) {
  assert(a_planar_separator_path.isIdSet());
  assert(this->isTagNew(a_planar_separator_path.getId()));
  id_mapping_m.insert(
      {a_planar_separator_path.getId(), std::move(a_planar_separator_path)});
}

inline void PlanarSeparatorPathGroup::addPlanarSeparatorPath(
    const PlanarSeparatorPath& a_planar_separator_path,
    const UnsignedIndex_t a_id) {
  assert(this->isTagNew(a_id));
  id_mapping_m.insert({a_id, a_planar_separator_path});
  id_mapping_m.at(a_id).setId(a_id);
}

inline void PlanarSeparatorPathGroup::addPlanarSeparatorPath(
    PlanarSeparatorPath&& a_planar_separator_path, const UnsignedIndex_t a_id) {
  assert(this->isTagNew(a_id));
  id_mapping_m.insert({a_id, std::move(a_planar_separator_path)});
  id_mapping_m.at(a_id).setId(a_id);
}

template <class StorageListType>
inline void PlanarSeparatorPathGroup::setPriorityOrder(
    const StorageListType& a_priority_order) {
  priority_order_by_id_m.resize(
      static_cast<std::size_t>(a_priority_order.size()));
  for (UnsignedIndex_t n = 0; n < priority_order_by_id_m.size(); ++n) {
    priority_order_by_id_m[n] = a_priority_order[n];
  }
  assert(allPrioritiesExist());
  this->setLinkingByPriorityList();
}
inline void PlanarSeparatorPathGroup::setPriorityOrder(
    const std::initializer_list<UnsignedIndex_t>& a_priority_order) {
  priority_order_by_id_m.resize(
      static_cast<std::size_t>(a_priority_order.size()));
  std::copy(a_priority_order.begin(), a_priority_order.end(),
            priority_order_by_id_m.begin());
  assert(allPrioritiesExist());
  this->setLinkingByPriorityList();
}
inline void PlanarSeparatorPathGroup::setPriorityOrder(
    const UnsignedIndex_t a_number_of_entries,
    const int* a_priority_order) {
  priority_order_by_id_m.resize(static_cast<std::size_t>(a_number_of_entries));
  for (UnsignedIndex_t n = 0; n < priority_order_by_id_m.size(); ++n) {
	assert(a_priority_order[n] >= 0);
    priority_order_by_id_m[n] = static_cast<UnsignedIndex_t>(a_priority_order[n]);
  }
  assert(allPrioritiesExist());
  this->setLinkingByPriorityList();
}

inline UnsignedIndex_t PlanarSeparatorPathGroup::getPriorityOrderSize(void) const{
  return static_cast<UnsignedIndex_t>(priority_order_by_id_m.size());
}
inline UnsignedIndex_t PlanarSeparatorPathGroup::getPriorityOrderTag(const UnsignedIndex_t a_index) const{
  assert(a_index < this->getPriorityOrderSize());
  return priority_order_by_id_m[a_index];
}

inline PlanarSeparatorPath&
PlanarSeparatorPathGroup::getReconstructionByPriority(
    const UnsignedIndex_t a_priority_index) {
  assert(a_priority_index < priority_order_by_id_m.size());
  assert(isTagKnown(priority_order_by_id_m[a_priority_index]));
  return id_mapping_m.at(priority_order_by_id_m[a_priority_index]);
}
inline PlanarSeparatorPath& PlanarSeparatorPathGroup::getReconstructionById(
    const UnsignedIndex_t a_id) {
  assert(isTagKnown(a_id));
  return id_mapping_m.at(a_id);
}

inline const PlanarSeparatorPath&
PlanarSeparatorPathGroup::getReconstructionByPriority(
    const UnsignedIndex_t a_priority_index) const {
  assert(a_priority_index < priority_order_by_id_m.size());
  assert(isTagKnown(priority_order_by_id_m[a_priority_index]));
  return id_mapping_m.at(priority_order_by_id_m[a_priority_index]);
}
inline const PlanarSeparatorPath&
PlanarSeparatorPathGroup::getReconstructionById(
    const UnsignedIndex_t a_id) const {
  assert(isTagKnown(a_id));
  return id_mapping_m.at(a_id);
}

inline PlanarSeparatorPath&
PlanarSeparatorPathGroup::getFirstReconstruction(void) {
  return this->getReconstructionByPriority(0);
}


inline const PlanarSeparatorPath&
PlanarSeparatorPathGroup::getFirstReconstruction(void) const {
  return this->getReconstructionByPriority(0);
}

inline const PlanarSeparatorPath&
PlanarSeparatorPathGroup::getCurrentReconstruction(void) const {
  return this->getReconstructionByPriority(0);
}

inline NullReconstruction PlanarSeparatorPathGroup::getNextReconstruction(
    void) const {
  return NullReconstruction();
}

inline void PlanarSeparatorPathGroup::setLinkingByPriorityList(void) {
  if (priority_order_by_id_m.empty()) {
    return;
  }

  PlanarSeparatorPath* current_reconstruction =
      &this->getReconstructionByPriority(0);
  for (UnsignedIndex_t n = 0; n < priority_order_by_id_m.size() - 1; ++n) {
    PlanarSeparatorPath* next_reconstruction =
        &this->getReconstructionByPriority(n + 1);
    current_reconstruction->setEdgeConnectivity(next_reconstruction);
    current_reconstruction = next_reconstruction;
  }
  current_reconstruction->setEdgeConnectivity(nullptr);
}

inline bool PlanarSeparatorPathGroup::isTagKnown(
    const UnsignedIndex_t a_tag) const {
  return id_mapping_m.find(a_tag) != id_mapping_m.end();
}
inline bool PlanarSeparatorPathGroup::isTagNew(
    const UnsignedIndex_t a_tag) const {
  return !this->isTagKnown(a_tag);
}

inline bool PlanarSeparatorPathGroup::allPrioritiesExist(void) const {
  for (const auto& index : priority_order_by_id_m) {
    if (isTagNew(index)) {
      // Prioritized something that is unknown.
      return false;
    }
  }
  return true;
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_GROUP_TPP_
