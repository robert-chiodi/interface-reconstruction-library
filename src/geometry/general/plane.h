// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PLANE_H_
#define SRC_GEOMETRY_GENERAL_PLANE_H_

#include <float.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <ostream>

#include "src/helpers/byte_buffer.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"
#include "src/helpers/serializer.h"
#include "src/parameters/constants.h"
#include "src/parameters/defined_types.h"

#include "src/geometry/general/normal.h"
#include "src/geometry/general/pt.h"

namespace IRL {

/// \brief Plane defined by \f$ \mathbf{n} \cdot \mathbf{x} - d = 0 \f$, where
/// the normal points from liquid to gas.
class Plane {
 public:
  /// \brief Default constructor.
  Plane(void);

  /// \brief Construct a plane given a Normal and a distance.
  Plane(const Normal& a_normal, const double a_distance);

  /// \brief Return a copy of the normal forming this plane.
  Normal& normal(void);

  /// \brief Return a const reference to the normal forming this plane.
  const Normal& normal(void) const;

  /// \brief Return value of the distance of the plane.
  double& distance(void);

  /// \brief Return const value to the distance of the plane.
  const double& distance(void) const;

  /// \brief Check if two planes are the same
  bool operator==(const Plane& a_other_plane) const;

  /// \brief Check if two planes are different
  bool operator!=(const Plane& a_other_plane) const;

  /// \brief Return signed distance from this plane to the supplied point, where
  /// negative is underneath the plane.
  template <class PtType>
  double signedDistanceToPoint(const PtType& a_pt) const;

  /// \brief Return size of the serialized point class in bytes.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize the point and store in the ByteBuffer
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the serialized normal and store.
  void unpackSerialized(ByteBuffer* a_buffer);

  Plane generateFlippedPlane(void) const;

  /// \brief Default destructor.
  ~Plane(void) = default;

 private:
  /// \brief Make sure normal have magnitude of 0.0 or 1.0.
  void checkValidNormal(void) const;

  Normal normal_m;    ///< \brief Normal for the plane.
  double distance_m;  ///< \brief Normal distance to the plane.
};

inline std::ostream& operator<<(std::ostream& out, const Plane& a_plane);

}  // namespace IRL

#include "src/geometry/general/plane.tpp"

#endif  // SRC_GEOMETRY_GENERAL_PLANE_H_
