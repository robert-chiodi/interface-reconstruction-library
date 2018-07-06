// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_NULL_RECONSTRUCTION_H_
#define SRC_PLANAR_RECONSTRUCTION_NULL_RECONSTRUCTION_H_

#include <iostream>

namespace IRL {

/// This is simply a blank class that is passed
/// when successive reconstructions are done
/// (such as in a LocalizedSeparator) to signify that
/// there are no more reconstructions to be done.
/// For instance, a PlanarSeparator does not have any
/// more reconstructions after it (its
/// getNextReconstruction(void) method should return
/// nothing). Instead, it returns this NullReconstruction,
/// which can have a function specialized for it where
/// the moments of the passed geometry object are simply
/// calculated and returned.
class NullReconstruction {
 public:
  NullReconstruction(void) = default;

  ~NullReconstruction(void) = default;

 private:
};

inline std::ostream& operator<<(std::ostream& out,
                                const NullReconstruction& a_reconstruction) {
  out << "This is a NullReconstruction" << std::endl;
  return out;
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_NULL_RECONSTRUCTION_H_
