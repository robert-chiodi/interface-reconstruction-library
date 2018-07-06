// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TO_GROUP_H_
#define SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TO_GROUP_H_

namespace IRL {

template <class CurrentReconstruction, class GroupReconstruction>
class JoinedReconstructionsToGroup {
 public:
  /// \brief Default constructor.
  JoinedReconstructionsToGroup(void);

  /// \brief Construct with providing a pointer to a localizer and separator.
  JoinedReconstructionsToGroup(
      const CurrentReconstruction* a_current_reconstruction_ptr,
      const GroupReconstruction* a_next_reconstruction_ptr);

  /// \brief Return the current reconstruction in JoinedReconstructions
  const CurrentReconstruction& getFirstReconstruction(void) const;

  /// \brief Return the next reconstruction in the JoinedReconstructions
  const GroupReconstruction& getSecondReconstruction(void) const;

  /// \brief Return the current reconstruction in JoinedReconstructions
  const auto& getCurrentReconstruction(void) const;

  /// \brief Return the next reconstruction in the JoinedReconstructions
  const auto& getNextReconstruction(void) const;

  /// \brief Default destructor;
  ~JoinedReconstructionsToGroup(void) = default;

 private:
  const CurrentReconstruction* current_reconstruction_m;
  const GroupReconstruction* next_reconstruction_m;
};

}  // namespace IRL

#include "src/planar_reconstruction/joined_reconstructions_to_group.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_TO_GROUP_H_
