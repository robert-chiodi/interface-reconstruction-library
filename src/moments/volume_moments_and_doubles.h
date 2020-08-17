// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_MOMENTS_AND_DOUBLES_H_
#define SRC_MOMENTS_VOLUME_MOMENTS_AND_DOUBLES_H_

#include "src/moments/volume_moments.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <UnsignedIndex_t kArrayLength>
class VolumeMomentsAndDoubles {
  using ArrayType = std::array<double, kArrayLength>;

 public:
  static constexpr UnsignedIndex_t data_length = kArrayLength;

  /// \brief Default constructor.
  VolumeMomentsAndDoubles(void);

  /// \brief Constructor that initializes volume and centroid.
  VolumeMomentsAndDoubles(const double a_volume, const Pt& a_centroid,
                          const ArrayType& a_initial_data);

  static VolumeMomentsAndDoubles fromRawDoublePointer(
      const double* a_list, const double* a_data_list);

  static VolumeMomentsAndDoubles fromScalarConstant(const double a_value);

  static VolumeMomentsAndDoubles fromScalarConstant(
      const double a_value, const double a_value_for_data);

  /// \brief Obtain un-normalized VolumeMoments from the supplied geometry.
  template <class GeometryType>
  static VolumeMomentsAndDoubles calculateMoments(GeometryType* a_geometry);

  /// \brief Return value of stored volume.
  Volume& volume(void);
  /// \brief Return const reference to stored volume.
  const Volume& volume(void) const;
  /// \brief Return copy of stored centroid.
  Pt& centroid(void);
  /// \brief Return const reference to stored centroid.
  const Pt& centroid(void) const;
  /// \brief Return reference to stored double array.
  ArrayType& data(void);
  /// \brief Return const reference to stored double array.
  const ArrayType& data(void) const;

  /// \brief Divide the centroid by the volume.
  void normalizeByVolume(void);

  /// \brief Multiply the centroid by the volume.
  void multiplyByVolume(void);

  /// \brief Overload += operator to update moments.
  VolumeMomentsAndDoubles& operator+=(const VolumeMomentsAndDoubles& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeMomentsAndDoubles& operator*=(const double a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeMomentsAndDoubles& operator/=(const double a_rhs);

  /// \brief Overload assignment to assign constant value to moments.
  VolumeMomentsAndDoubles& operator=(const double a_value);

  /// \brief Default destructor.
  ~VolumeMomentsAndDoubles(void) = default;

 private:
  /// \brief Construct VolumeMoments from a list of doubles.
  ///
  /// Construct VolumeMoments from a list of 4 doubles. The necessary order
  /// is volume, centroid_x, centroid_y, centroid_z.
  explicit VolumeMomentsAndDoubles(const double* a_list,
                                   const double* a_data_list);

  /// \brief Construct that initializes volume/centroid as a value.
  explicit VolumeMomentsAndDoubles(const double a_value,
                                   const double a_value_for_data);

  VolumeMoments volume_moments_m;
  std::array<double, kArrayLength> data_moments_m;
};

template <UnsignedIndex_t kArrayLength>
std::ostream& operator<<(
    std::ostream& out,
    const VolumeMomentsAndDoubles<kArrayLength>& a_volume_moments_and_doubles);

template <UnsignedIndex_t kArrayLength>
VolumeMomentsAndDoubles<kArrayLength> operator+(
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad1,
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad2);

template <UnsignedIndex_t kArrayLength>
VolumeMomentsAndDoubles<kArrayLength> operator-(
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad1,
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad2);

template <UnsignedIndex_t kArrayLength>
VolumeMomentsAndDoubles<kArrayLength> operator*(
    const double a_multiplier,
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad);

template <UnsignedIndex_t kArrayLength>
VolumeMomentsAndDoubles<kArrayLength> operator*(
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad,
    const double a_multiplier);

}  // namespace IRL

#include "src/moments/volume_moments_and_doubles.tpp"

#endif  // SRC_MOMENTS_VOLUME_MOMENTS_AND_DOUBLES_H_
