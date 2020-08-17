// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_ACCUMULATE_WRAPPER_H_
#define SRC_MOMENTS_ACCUMULATE_WRAPPER_H_

namespace IRL {

/// \brief This is essentially a class that delays the evaluation of Tagged
/// moments during generic_cutting with a ReconstructionLink for an eventual
/// collection. As an example, this class needs to be used when cutting by a
/// LocalizedSeparatorPath for tagged moments, since the first cutting operation
/// (by the Localizer), is not what will dictate the tagging of the moments, but
/// instead the second (by the PlanarSeparatorPath) does. The use of the tagged
/// moments that are going to be returned need to be delayed and accumulated. We
/// would then return them as auto moments =
/// getVolumeMoments<AccumulateWrapper<TaggedAccumulatedVolumeMoments<VolumeMoments>>>
/// (a_volume_to_cut, a_LocalizedSeparatorPath);
template <class WrappedType>
class AccumulateWrapper : public WrappedType {
 public:
  using WrappedType::WrappedType;

  using contained_type = WrappedType;

 private:
};
}  // namespace IRL

#endif  // SRC_MOMENTS_ACCUMULATE_WRAPPER_H_
