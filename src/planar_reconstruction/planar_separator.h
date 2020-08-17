// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_H_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_H_

#include <ostream>

#include "src/data_structures/small_vector.h"
#include "src/geometry/general/plane.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/serializer.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/null_reconstruction.h"
#include "src/planar_reconstruction/planar_reconstruction.h"

namespace IRL {

class PlanarSeparator {
  using PlanarReconstructionBase =
      PlanarReconstruction<global_constants::MAX_PLANAR_LOCALIZER_PLANES>;
  using iterator = PlanarReconstructionBase::iterator;
  using const_iterator = PlanarReconstructionBase::const_iterator;
  friend std::ostream& operator<<(std::ostream& out,
                                  const PlanarSeparator& a_reconstruction);

 public:
  /// \brief Default constructor, initialize one plane that is far below
  /// everything.
  PlanarSeparator(void);

  static PlanarSeparator fromOnePlane(const Plane& a_plane);

  static PlanarSeparator fromTwoPlanes(const Plane& a_plane_0,
                                       const Plane& a_plane_1,
                                       const double a_flip_indicator);

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

  /// \brief Add a new plane to the reconstruction
  void addPlane(const Plane& a_plane);

  /// \brief Set distances for planes by given a pointer to
  /// an array of doubles this->getNumberOfPlanes() long.
  template <class ArrayType>
  void setDistances(const ArrayType& a_distances);

  /// \brief Remove the plane given by index `a_p`
  void removePlane(const UnsignedIndex_t a_p);

  /// \brief Return value of `flip_m`
  double flip(void) const;

  /// \brief Set flip directly
  void setFlip(const double a_flip_value);

  /// \brief Make it so regular cutting is used
  /// (liquid phase, below plane, found by cutting).
  void doNotFlipCutting(void);

  /// \brief Make it so that cutting by this reconstruction will be flipped
  /// (gas phase found through cutting).
  void flipCutting(void);

  /// \brief Return if cutting for gas phase is needed.
  bool isFlipped(void) const;

  /// \brief Return if cutting for gas phase is needed.
  bool isNotFlipped(void) const;

  /// \brief Set reconstruction to cause single-phase cell.
  void zeroPlanes(void);

  PlanarSeparator& getCurrentReconstruction(void);
  const PlanarSeparator& getCurrentReconstruction(void) const;
  static constexpr NullReconstruction getNextReconstruction(void);

  /// \brief Return size of the serialized PlanarSeparator.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the planes and flip_cut_m.
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the planes and flip_cut_m.
  void unpackSerialized(ByteBuffer* a_buffer);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor
  ~PlanarSeparator(void) = default;

 private:
  /// \brief Constructor when given a single plane
  explicit PlanarSeparator(const Plane& a_plane);

  /// \brief Constructor when given two planes
  PlanarSeparator(const Plane& a_plane_0, const Plane& a_plane_1,
                  const double a_flip_indicator);

  PlanarReconstructionBase reconstruction_m;
  double flip_cut_m;  ///< \brief Used to flip the phase cut for.
};

inline std::ostream& operator<<(std::ostream& out,
                                const PlanarSeparator& a_reconstruction);

}  // namespace IRL

#include "src/planar_reconstruction/planar_separator.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_H_
