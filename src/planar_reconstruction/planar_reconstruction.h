// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_H_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_H_

#include <algorithm>
#include <ostream>
#include <vector>

#include "src/data_structures/small_vector.h"
#include "src/geometry/general/plane.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/serializer.h"

namespace IRL {

template <UnsignedIndex_t kStackPlanes>
class PlanarReconstruction {
  friend class PlanarSeparator;
  friend class PlanarLocalizer;
  using iterator = typename SmallVector<Plane, kStackPlanes>::iterator;
  using const_iterator =
      typename SmallVector<Plane, kStackPlanes>::const_iterator;

 public:
  /// \brief Default constructor.
  PlanarReconstruction(void) = default;

  static PlanarReconstruction fromOnePlane(const Plane& a_plane);

  static PlanarReconstruction fromTwoPlanes(const Plane& a_plane_0,
                                            const Plane& a_plane_1);

  /// \brief Return the number of planes
  /// used for the reconstruction.
  UnsignedIndex_t getNumberOfPlanes(void) const;

  /// \brief Reserve space for a number of planes that will (presumably) be
  /// added later.
  void reserve(const UnsignedIndex_t a_number_of_future_planes);

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

  /// \brief Set distances to the planes in the reconstruction
  template <class ArrayType>
  void setDistances(const ArrayType& a_distances);

  /// \brief Remove the plane given by index `a_p`
  void removePlane(const UnsignedIndex_t a_p);

  /// \brief Return size of the serialized PlanarReconstruction.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the planes.
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the planes and store.
  void unpackSerialized(ByteBuffer* a_buffer);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~PlanarReconstruction(void) = default;

 private:
  /// \brief Construct with one given plane
  explicit PlanarReconstruction(const Plane& a_plane);

  /// \brief Construct with two given planes
  PlanarReconstruction(const Plane& a_plane_0, const Plane& a_plane_1);

  void checkIfStaticAllocationExceeded(void) const;

  /// \brief Planes in reconstruction.
  SmallVector<Plane, kStackPlanes> planes_m;
};

template <UnsignedIndex_t kStackPlanes>
inline std::ostream& operator<<(
    std::ostream& out,
    const PlanarReconstruction<kStackPlanes>& a_reconstruction);

}  // namespace IRL

#include "src/planar_reconstruction/planar_reconstruction.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_H_
