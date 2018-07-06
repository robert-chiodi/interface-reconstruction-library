// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_ACCUMULATED_LISTED_VOLUME_MOMENTS_H_
#define SRC_MOMENTS_ACCUMULATED_LISTED_VOLUME_MOMENTS_H_

#include "src/moments/accumulated_volume_moments.h"
#include "src/moments/listed_volume_moments.h"

namespace IRL {

template <class VolumeMomentsType>
class AccumulatedListedVolumeMoments
    : public AccumulatedVolumeMoments<ListedVolumeMoments<VolumeMomentsType>> {
 public:
  using contained_type = VolumeMomentsType;
  /// \brief Default constructor.
  AccumulatedListedVolumeMoments(void) = default;

  /// \brief Default destructor.
  ~AccumulatedListedVolumeMoments(void) = default;

 private:
};
}  // namespace IRL

#endif  // SRC_MOMENTS_ACCUMULATED_LISTED_VOLUME_MOMENTS_H_
