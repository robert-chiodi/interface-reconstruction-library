// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_LOCALIZER_H_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_LOCALIZER_H_

#include <ostream>

#include "src/geometry/general/plane.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/serializer.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/null_reconstruction.h"
#include "src/planar_reconstruction/planar_reconstruction.h"

namespace IRL {

/// \brief A planar representation of a convex polyhedron to localize
/// integrations.
class PlanarLocalizer {
  using PlanarReconstructionBase =
      PlanarReconstruction<global_constants::MAX_PLANAR_LOCALIZER_PLANES>;
  using iterator = PlanarReconstructionBase::iterator;
  using const_iterator = PlanarReconstructionBase::const_iterator;
  friend std::ostream& operator<<(std::ostream& out,
                                  const PlanarLocalizer& a_reconstruction);

 public:
  /// \brief Default constructor
  PlanarLocalizer(void) = default;

  static PlanarLocalizer fromOnePlane(const Plane& a_plane);

  static PlanarLocalizer fromTwoPlanes(const Plane& a_plane_0,
                                       const Plane& a_plane_1);

  PlanarLocalizer(const PlanarLocalizer& a_other) = default;
  PlanarLocalizer& operator=(const PlanarLocalizer& a_other) = default;

  /// \brief Return the number of planes
  /// used for the reconstruction.
  UnsignedIndex_t getNumberOfPlanes(void) const;

  /// \brief Shrink/enlarge vector of planes to given number.
  void setNumberOfPlanes(const UnsignedIndex_t a_number_of_future_planes);

  /// \brief Directly set the number of planes to 0.
  void zeroNumberOfPlanes(void);

  /// \brief Overload `[]` to access planes_m through reference.
  Plane& operator[](const UnsignedIndex_t a_p);

  /// \brief Overload `[]` to access planes_m through const reference.
  const Plane& operator[](const UnsignedIndex_t a_p) const;

  /// \brief Add a new plane to the list
  void addPlane(const Plane& a_plane);

  /// \brief Insert new plane at start of planes.
  void addBeginningPlane(const Plane& a_plane);

  /// \brief Remove the plane given by index `a_p`
  void removePlane(const UnsignedIndex_t a_p);

  /// \brief Since localizers are always convex, never flip.
  static constexpr double flip(void);

  /// \brief Since localizers are always convex, never flipped.
  static constexpr bool isFlipped(void);

  /// \brief Return if cutting for gas phase is needed.
  static constexpr bool isNotFlipped(void);

  const PlanarLocalizer& getCurrentReconstruction(void) const;
  static constexpr NullReconstruction getNextReconstruction(void);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Return size of the serialized PlanarLocalizer.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the planes.
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the planes.
  void unpackSerialized(ByteBuffer* a_buffer);

  /// \brief Default destructor
  ~PlanarLocalizer(void) = default;

 private:
  /// \brief Constructor when given a single plane
  explicit PlanarLocalizer(const Plane& a_plane);

  /// \brief Constructor when given two planes
  PlanarLocalizer(const Plane& a_plane_0, const Plane& a_plane_1);

  PlanarReconstructionBase reconstruction_m;
};

inline std::ostream& operator<<(std::ostream& out,
                                const PlanarLocalizer& a_reconstruction);

}  // namespace IRL

#include "src/planar_reconstruction/planar_localizer.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_LOCALIZER_H_
