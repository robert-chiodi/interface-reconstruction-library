// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_GROUP_H_
#define IRL_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_GROUP_H_

#include <array>
#include <initializer_list>
#include <vector>

#include "irl/data_structures/unordered_map.h"
#include "irl/parameters/defined_types.h"
#include "irl/planar_reconstruction/null_reconstruction.h"
#include "irl/planar_reconstruction/planar_separator_path.h"

namespace IRL {

/// \brief This is a class that helps organize PlanarSeparatorPath
/// objects and facilitates their control.
/// It will continue to store all PlanarSeparatorPath objects that were added to
/// it, but in most routines only the active ones will be used. Here, active
/// means those that are present in priority_order_by_id_m.
class PlanarSeparatorPathGroup {
 public:
  PlanarSeparatorPathGroup(void) = default;

  void addPlanarSeparatorPath(
      const PlanarSeparatorPath& a_planar_separator_path);
  void addPlanarSeparatorPath(PlanarSeparatorPath&& a_planar_separator_path);

  void addPlanarSeparatorPath(
      const PlanarSeparatorPath& a_planar_separator_path,
      const UnsignedIndex_t a_id);
  void addPlanarSeparatorPath(PlanarSeparatorPath&& a_planar_separator_path,
                              const UnsignedIndex_t a_id);

  /// \brief StorageListType must have a []() operator and size() method.
  template <class StorageListType>
  void setPriorityOrder(const StorageListType& a_priority_order);
  void setPriorityOrder(
      const std::initializer_list<UnsignedIndex_t>& a_priority_order);
  void setPriorityOrder(const UnsignedIndex_t a_number_of_entries,
                        const int* a_priority_order);

  UnsignedIndex_t getPriorityOrderSize(void) const;
  UnsignedIndex_t getPriorityOrderTag(const UnsignedIndex_t a_index) const;

  PlanarSeparatorPath& getReconstructionByPriority(
      const UnsignedIndex_t a_priority_index);
  const PlanarSeparatorPath& getReconstructionByPriority(
      const UnsignedIndex_t a_priority_index) const;
  const PlanarSeparatorPath& getReconstructionById(
      const UnsignedIndex_t a_id) const;
  PlanarSeparatorPath& getReconstructionById(const UnsignedIndex_t a_id);

  PlanarSeparatorPath& getFirstReconstruction(void);
  const PlanarSeparatorPath& getFirstReconstruction(void) const;

  const PlanarSeparatorPath& getCurrentReconstruction(void) const;
  NullReconstruction getNextReconstruction(void) const;

  ~PlanarSeparatorPathGroup(void) = default;

 private:
  void setLinkingByPriorityList(void);

  bool isTagKnown(const UnsignedIndex_t a_tag) const;
  bool isTagNew(const UnsignedIndex_t a_tag) const;
  bool allPrioritiesExist(void) const;

  std::vector<UnsignedIndex_t> priority_order_by_id_m;
  IRL::unordered_map<UnsignedIndex_t, PlanarSeparatorPath> id_mapping_m;
};

}  // namespace IRL

#include "irl/planar_reconstruction/planar_separator_path_group.tpp"

#endif // IRL_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_GROUP_H_
