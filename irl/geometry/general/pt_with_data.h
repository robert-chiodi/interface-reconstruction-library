// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_PT_WITH_DATA_H_
#define IRL_GEOMETRY_GENERAL_PT_WITH_DATA_H_

#include <array>

#include "irl/geometry/general/pt.h"
#include "irl/helpers/helper.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class Derived, class AttachedDataType, class ScalarType>
class PtWithDataCommon {
  Derived& getDerived(void);
  const Derived& getDerived(void) const;

 public:
  using contained_type = AttachedDataType;
  using value_type = ScalarType;

  PtWithDataCommon(void) = default;

  PtWithDataCommon(const PtBase<ScalarType>& a_pt,
                   const AttachedDataType& a_data);

  explicit PtWithDataCommon(const PtBase<ScalarType>& a_pt);

  /// \brief Provides access to underlying base_pt location.
  ScalarType& operator[](const UnsignedIndex_t a_index);
  const ScalarType& operator[](const UnsignedIndex_t a_index) const;

  PtBase<ScalarType>& getPt(void);
  const PtBase<ScalarType>& getPt(void) const;

  AttachedDataType& getData(void);
  const AttachedDataType& getData(void) const;

  ~PtWithDataCommon(void) = default;

 private:
  PtBase<ScalarType> base_pt_m;
  AttachedDataType data_m;
};

template <class GradientDataType, class ScalarType>
class PtWithGradientBase
    : public PtWithDataCommon<PtWithGradientBase<GradientDataType, ScalarType>,
                              std::array<GradientDataType, 3>, ScalarType> {
 public:
  using gradient_type = GradientDataType;

  PtWithGradientBase(void);

  using PtWithDataCommon<PtWithGradientBase<GradientDataType, ScalarType>,
                         std::array<GradientDataType, 3>,
                         ScalarType>::PtWithDataCommon;

  explicit PtWithGradientBase(const PtBase<ScalarType>& a_pt);

  PtWithGradientBase& operator-(void);
  PtWithGradientBase& operator=(const PtBase<ScalarType>& a_pt);
  PtWithGradientBase& operator=(const PtWithGradientBase& a_pt);
  PtWithGradientBase& operator+=(const PtWithGradientBase& a_other_pt);

  ~PtWithGradientBase(void) = default;
};

template <class GradientDataType, class ScalarType>
const PtWithGradientBase<GradientDataType, ScalarType> operator*(
    const ScalarType a_rhs,
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt);
template <class GradientDataType, class ScalarType>
const PtWithGradientBase<GradientDataType, ScalarType> operator*(
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt,
    const ScalarType a_rhs);
template <class GradientDataType, class ScalarType>
const PtWithGradientBase<GradientDataType, ScalarType> operator/(
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt,
    const ScalarType a_rhs);
template <class GradientDataType, class ScalarType>
const PtWithGradientBase<GradientDataType, ScalarType> operator+(
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt1,
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt2);
template <class GradientDataType, class ScalarType>
const PtWithGradientBase<GradientDataType, ScalarType> operator-(
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt1,
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt2);
template <class GradientDataType, class ScalarType>
const PtWithGradientBase<GradientDataType, ScalarType> operator-(
    const PtWithGradientBase<GradientDataType, ScalarType>& a_pt);

template <class GradientDataType>
using PtWithGradient = PtWithGradientBase<GradientDataType, double>;

template <class FunctorType, UnsignedIndex_t kArrayLength>
class PtWithDoublesStatelessFunctor
    : public PtWithDataCommon<
          PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
          std::array<double, kArrayLength>, double> {
  using ArrayType = std::array<double, kArrayLength>;

 public:
  static constexpr UnsignedIndex_t data_length = kArrayLength;

  PtWithDoublesStatelessFunctor(void) = default;

  using PtWithDataCommon<
      PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
      std::array<double, kArrayLength>, double>::PtWithDataCommon;

  explicit PtWithDoublesStatelessFunctor(const Pt& a_pt);
  PtWithDoublesStatelessFunctor& operator=(const Pt& a_pt);

  static PtWithDoublesStatelessFunctor fromEdgeIntersection(
      const PtWithDoublesStatelessFunctor& a_pt_0, const double a_dist_0,
      const PtWithDoublesStatelessFunctor& a_pt_1, const double a_dist_1);

  PtWithDoublesStatelessFunctor& operator+=(
      const PtWithDoublesStatelessFunctor& a_other_pt);

  PtWithDoublesStatelessFunctor& operator/=(const double a_double);

  ~PtWithDoublesStatelessFunctor(void) = default;
};

template <class FunctorType, UnsignedIndex_t kArrayLength>
std::ostream& operator<<(
    std::ostream& out,
    const PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>&
        a_pt_with_data);

}  // namespace IRL

#include "irl/geometry/general/pt_with_data.tpp"

#endif  // IRL_GEOMETRY_GENERAL_PT_WITH_DATA_H_
