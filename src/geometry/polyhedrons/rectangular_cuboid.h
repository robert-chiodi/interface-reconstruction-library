// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_RECTANGULAR_CUBOID_H_
#define SRC_GEOMETRY_POLYHEDRONS_RECTANGULAR_CUBOID_H_

#include <float.h>

#include <cmath>

#include "src/geometry/general/normal.h"
#include "src/geometry/general/pt.h"
#include "src/geometry/general/pt_list.h"
#include "src/geometry/general/stored_vertex_access.h"
#include "src/geometry/polyhedrons/hexahedron.h"
#include "src/helpers/mymath.h"
#include "src/moments/volume_moments.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/planar_localizer.h"

namespace IRL {

/// \brief A rectangular cuboid.
template <class Derived, class VertexType>
class RectangularCuboidSpecialization
    : public HexahedronSpecialization<Derived, VertexType> {
 public:
  /// \brief Calculate and return side length for `a_dimension`.
  double calculateSideLength(const UnsignedIndex_t a_dimension) const;

  /// \brief Returns a planar reconstruction that is equivalent to the
  /// RectangularCuboidCommon.
  PlanarLocalizer getLocalizer(void) const;

  /// \brief Calculate and return the volume of the rectangular cuboid.
  inline double calculateVolume(void) const;

  /// \brief Calculate and return the centroid of the rectangular cuboid.
  inline Pt calculateCentroid(void) const;

  /// \brief Calculate and return volume weighted VolumeMoments.
  inline VolumeMoments calculateMoments() const;

 private:
};

template <class VertexType>
class StoredRectangularCuboid
    : public StoredVertexAccess<StoredRectangularCuboid<VertexType>, VertexType,
                                8>,
      public RectangularCuboidSpecialization<
          StoredRectangularCuboid<VertexType>, VertexType> {
  friend StoredVertexAccess<StoredRectangularCuboid<VertexType>, VertexType, 8>;

 public:
  using StoredVertexAccess<StoredRectangularCuboid<VertexType>, VertexType,
                           8>::StoredVertexAccess;

  static StoredRectangularCuboid fromBoundingPts(const Pt& a_min_point,
                                                 const Pt& a_max_point);

  StoredRectangularCuboid(void) = default;

 private:
  StoredRectangularCuboid(const Pt& a_min_point, const Pt& a_max_point);
};

// Predefined types
using RectangularCuboid = StoredRectangularCuboid<Pt>;

extern const RectangularCuboid unit_cell;

}  // namespace IRL

#include "src/geometry/polyhedrons/rectangular_cuboid.tpp"

#endif  // SRC_GEOMETRY_POLYHEDRONS_RECTANGULAR_CUBOID_H_
