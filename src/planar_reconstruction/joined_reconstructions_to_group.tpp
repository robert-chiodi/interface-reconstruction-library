// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TO_GROUP_TPP_
#define SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TO_GROUP_TPP_

namespace IRL {

template <class CurrentReconstruction, class GroupReconstruction>
inline JoinedReconstructionsToGroup<
    CurrentReconstruction,
    GroupReconstruction>::JoinedReconstructionsToGroup(void)
    : current_reconstruction_m(nullptr), next_reconstruction_m(nullptr) {}

template <class CurrentReconstruction, class GroupReconstruction>
inline JoinedReconstructionsToGroup<CurrentReconstruction,
                                    GroupReconstruction>::
    JoinedReconstructionsToGroup(
        const CurrentReconstruction* a_current_reconstruction_ptr,
        const GroupReconstruction* a_next_reconstruction_ptr)
    : current_reconstruction_m(a_current_reconstruction_ptr),
      next_reconstruction_m(a_next_reconstruction_ptr) {
  assert(a_current_reconstruction_ptr != nullptr);
  assert(a_next_reconstruction_ptr != nullptr);
}

template <class CurrentReconstruction, class GroupReconstruction>
inline const CurrentReconstruction& JoinedReconstructionsToGroup<
    CurrentReconstruction, GroupReconstruction>::getFirstReconstruction(void)
    const {
  assert(current_reconstruction_m != nullptr);
  return *current_reconstruction_m;
}

template <class CurrentReconstruction, class GroupReconstruction>
inline const GroupReconstruction& JoinedReconstructionsToGroup<
    CurrentReconstruction, GroupReconstruction>::getSecondReconstruction(void)
    const {
  assert(next_reconstruction_m != nullptr);
  return *next_reconstruction_m;
}

template <class CurrentReconstruction, class GroupReconstruction>
inline const auto& JoinedReconstructionsToGroup<
    CurrentReconstruction, GroupReconstruction>::getCurrentReconstruction(void)
    const {
  assert(current_reconstruction_m != nullptr);
  return current_reconstruction_m->getCurrentReconstruction();
}

template <class CurrentReconstruction, class GroupReconstruction>
inline const auto& JoinedReconstructionsToGroup<
    CurrentReconstruction, GroupReconstruction>::getNextReconstruction(void)
    const {
  assert(next_reconstruction_m != nullptr);
  return next_reconstruction_m->getFirstReconstruction();
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TO_GROUP_TPP_
