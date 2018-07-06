// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PT_WITH_DATA_H_
#define SRC_GEOMETRY_GENERAL_PT_WITH_DATA_H_

#include <array>

#include "src/geometry/general/pt.h"
#include "src/helpers/helper.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class Derived, class AttachedDataType>
class PtWithDataCommon {
  Derived& getDerived(void);
  const Derived& getDerived(void) const;

 public:
  using contained_type = AttachedDataType;

  PtWithDataCommon(void) = default;

  PtWithDataCommon(const Pt& a_pt, const AttachedDataType& a_data);

  explicit PtWithDataCommon(const Pt& a_pt);

  /// \brief Provides access to underlying base_pt location.
  double& operator[](const UnsignedIndex_t a_index);
  const double& operator[](const UnsignedIndex_t a_index) const;

  Pt& getPt(void);
  const Pt& getPt(void) const;

  AttachedDataType& getData(void);
  const AttachedDataType& getData(void) const;

  ~PtWithDataCommon(void) = default;

 private:
  Pt base_pt_m;
  AttachedDataType data_m;
};

template <class FunctorType, UnsignedIndex_t kArrayLength>
class PtWithDoublesStatelessFunctor
    : public PtWithDataCommon<
          PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
          std::array<double, kArrayLength>> {
  using ArrayType = std::array<double, kArrayLength>;

 public:
  static constexpr UnsignedIndex_t data_length = kArrayLength;

  PtWithDoublesStatelessFunctor(void) = default;

  using PtWithDataCommon<
      PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
      std::array<double, kArrayLength>>::PtWithDataCommon;

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

#include "src/geometry/general/pt_with_data.tpp"

#endif  // SRC_GEOMETRY_GENERAL_PT_WITH_DATA_H_
