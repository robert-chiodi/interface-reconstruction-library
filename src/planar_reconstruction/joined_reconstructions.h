// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_H_
#define SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_H_

#include "src/planar_reconstruction/null_reconstruction.h"

namespace IRL {

/// \brief Class that ties together two reconstructions. Will first
/// cut by and use the current reconstruction, and then pass
/// onwards to the next reconstruction.
template <class CurrentReconstruction, class NextReconstruction>
class JoinedReconstructions {
 public:
  /// \brief Default constructor.
  JoinedReconstructions(void);

  /// \brief Construct with providing a pointer to a localizer and separator.
  JoinedReconstructions(
      const CurrentReconstruction* a_current_reconstruction_ptr,
      const NextReconstruction* a_next_reconstruction_ptr);

  /// \brief Return the current reconstruction in JoinedReconstructions
  const CurrentReconstruction& getFirstReconstruction(void) const;

  /// \brief Return the next reconstruction in the JoinedReconstructions
  const NextReconstruction& getSecondReconstruction(void) const;

  /// \brief Return the current reconstruction in JoinedReconstructions
  const auto& getCurrentReconstruction(void) const;

  /// \brief Return the next reconstruction in the JoinedReconstructions
  const NextReconstruction& getNextReconstruction(void) const;

  /// \brief Default destructor;
  ~JoinedReconstructions(void) = default;

 private:
  const CurrentReconstruction* current_reconstruction_m;
  const NextReconstruction* next_reconstruction_m;
};

template <class CurrentReconstruction>
class JoinedReconstructions<CurrentReconstruction, NullReconstruction> {
  using NextReconstruction = NullReconstruction;

 public:
  /// \brief Default constructor.
  JoinedReconstructions(void);

  /// \brief Construct with providing a pointer to a localizer and separator.
  JoinedReconstructions(
      const CurrentReconstruction* a_current_reconstruction_ptr);

  /// \brief Return the first reconstruction (CurrentReconstruction)
  const CurrentReconstruction& getFirstReconstruction(void) const;

  /// \brief Return the second reconstruction (NextReconstruction)
  NullReconstruction getSecondReconstruction(void) const;

  /// \brief Return the current reconstruction of the first reconstruction
  const auto& getCurrentReconstruction(void) const;

  /// \brief Return the second reconstruction.
  NullReconstruction getNextReconstruction(void) const;

  /// \brief Default destructor;
  ~JoinedReconstructions(void) = default;

 private:
  const CurrentReconstruction* current_reconstruction_m;
};

}  // namespace IRL

#include "src/planar_reconstruction/joined_reconstructions.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_JOINED_RECONSTRUCTIONS_H_
