// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TPP_
#define SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TPP_

namespace IRL {

template <class CurrentReconstruction, class NextReconstruction>
inline JoinedReconstructions<CurrentReconstruction,
                             NextReconstruction>::JoinedReconstructions(void)
    : current_reconstruction_m(nullptr), next_reconstruction_m(nullptr) {}

template <class CurrentReconstruction, class NextReconstruction>
inline JoinedReconstructions<CurrentReconstruction, NextReconstruction>::
    JoinedReconstructions(
        const CurrentReconstruction* a_current_reconstruction_ptr,
        const NextReconstruction* a_next_reconstruction_ptr)
    : current_reconstruction_m(a_current_reconstruction_ptr),
      next_reconstruction_m(a_next_reconstruction_ptr) {
  assert(a_current_reconstruction_ptr != nullptr);
  assert(a_next_reconstruction_ptr != nullptr);
}

template <class CurrentReconstruction, class NextReconstruction>
inline const CurrentReconstruction&
JoinedReconstructions<CurrentReconstruction,
                      NextReconstruction>::getFirstReconstruction(void) const {
  assert(current_reconstruction_m != nullptr);
  return *current_reconstruction_m;
}

template <class CurrentReconstruction, class NextReconstruction>
inline const NextReconstruction&
JoinedReconstructions<CurrentReconstruction,
                      NextReconstruction>::getSecondReconstruction(void) const {
  assert(next_reconstruction_m != nullptr);
  return *next_reconstruction_m;
}

template <class CurrentReconstruction, class NextReconstruction>
inline const auto& JoinedReconstructions<
    CurrentReconstruction, NextReconstruction>::getCurrentReconstruction(void)
    const {
  assert(current_reconstruction_m != nullptr);
  return current_reconstruction_m->getCurrentReconstruction();
}

template <class CurrentReconstruction, class NextReconstruction>
inline const NextReconstruction&
JoinedReconstructions<CurrentReconstruction,
                      NextReconstruction>::getNextReconstruction(void) const {
  assert(next_reconstruction_m != nullptr);
  return *next_reconstruction_m;
}

//******************************************************************* //
//     Null Reconstruction specialization below this
//******************************************************************* //
template <class CurrentReconstruction>
inline JoinedReconstructions<CurrentReconstruction,
                             NullReconstruction>::JoinedReconstructions(void)
    : current_reconstruction_m(nullptr) {}

template <class CurrentReconstruction>
inline JoinedReconstructions<CurrentReconstruction, NullReconstruction>::
    JoinedReconstructions(
        const CurrentReconstruction* a_current_reconstruction_ptr)
    : current_reconstruction_m(a_current_reconstruction_ptr) {
  assert(a_current_reconstruction_ptr != nullptr);
}

template <class CurrentReconstruction>
inline const CurrentReconstruction&
JoinedReconstructions<CurrentReconstruction,
                      NullReconstruction>::getFirstReconstruction(void) const {
  assert(current_reconstruction_m != nullptr);
  return *current_reconstruction_m;
}

template <class CurrentReconstruction>
inline NullReconstruction
JoinedReconstructions<CurrentReconstruction,
                      NullReconstruction>::getSecondReconstruction(void) const {
  return NullReconstruction();
}

template <class CurrentReconstruction>
inline const auto& JoinedReconstructions<
    CurrentReconstruction, NullReconstruction>::getCurrentReconstruction(void)
    const {
  assert(current_reconstruction_m != nullptr);
  return current_reconstruction_m->getCurrentReconstruction();
}

template <class CurrentReconstruction>
inline NullReconstruction
JoinedReconstructions<CurrentReconstruction,
                      NullReconstruction>::getNextReconstruction(void) const {
  return NullReconstruction();
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TPP_
