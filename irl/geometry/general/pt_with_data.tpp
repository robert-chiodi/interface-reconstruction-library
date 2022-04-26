// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_PT_WITH_DATA_TPP_
#define IRL_GEOMETRY_GENERAL_PT_WITH_DATA_TPP_

namespace IRL {

template <class Derived, class AttachedDataType>
Derived& PtWithDataCommon<Derived, AttachedDataType>::getDerived(void) {
  return static_cast<Derived&>(*this);
}
template <class Derived, class AttachedDataType>
const Derived& PtWithDataCommon<Derived, AttachedDataType>::getDerived(
    void) const {
  return static_cast<const Derived&>(*this);
}

template <class Derived, class AttachedDataType>
PtWithDataCommon<Derived, AttachedDataType>::PtWithDataCommon(
    const Pt& a_pt, const AttachedDataType& a_data)
    : base_pt_m(a_pt), data_m(a_data) {}

template <class Derived, class AttachedDataType>
PtWithDataCommon<Derived, AttachedDataType>::PtWithDataCommon(const Pt& a_pt)
    : base_pt_m(a_pt) {}

template <class Derived, class AttachedDataType>
double& PtWithDataCommon<Derived, AttachedDataType>::operator[](
    const UnsignedIndex_t a_index) {
  return base_pt_m[a_index];
}
template <class Derived, class AttachedDataType>
const double& PtWithDataCommon<Derived, AttachedDataType>::operator[](
    const UnsignedIndex_t a_index) const {
  return base_pt_m[a_index];
}

template <class Derived, class AttachedDataType>
Pt& PtWithDataCommon<Derived, AttachedDataType>::getPt(void) {
  return base_pt_m;
}

template <class Derived, class AttachedDataType>
const Pt& PtWithDataCommon<Derived, AttachedDataType>::getPt(void) const {
  return base_pt_m;
}

template <class Derived, class AttachedDataType>
AttachedDataType& PtWithDataCommon<Derived, AttachedDataType>::getData(void) {
  return data_m;
}
template <class Derived, class AttachedDataType>
const AttachedDataType& PtWithDataCommon<Derived, AttachedDataType>::getData(
    void) const {
  return data_m;
}

template <class FunctorType, UnsignedIndex_t kArrayLength>
PtWithDoublesStatelessFunctor<
    FunctorType, kArrayLength>::PtWithDoublesStatelessFunctor(const Pt& a_pt)
    : PtWithDataCommon<PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
                       std::array<double, kArrayLength>>(a_pt) {
  this->getData().fill(0.0);
}

template <class FunctorType, UnsignedIndex_t kArrayLength>
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>&
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>::operator=(
    const Pt& a_pt) {
  this->getPt() = a_pt;
  this->getData().fill(0.0);
  return (*this);
}

template <class FunctorType, UnsignedIndex_t kArrayLength>
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>::fromEdgeIntersection(
    const PtWithDoublesStatelessFunctor& a_pt_0, const double a_dist_0,
    const PtWithDoublesStatelessFunctor& a_pt_1, const double a_dist_1) {
  FunctorType func;
  return func(a_pt_0, a_dist_0, a_pt_1, a_dist_1);
}

template <class FunctorType, UnsignedIndex_t kArrayLength>
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>&
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>::operator+=(
    const PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>&
        a_other_pt) {
  this->getPt() += a_other_pt.getPt();
  const auto& other_data = a_other_pt.getData();
  for (UnsignedIndex_t n = 0; n < other_data.size(); ++n) {
    this->getData()[n] += other_data[n];
  }
  return (*this);
}

template <class FunctorType, UnsignedIndex_t kArrayLength>
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>&
PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>::operator/=(
    const double a_double) {
  this->getPt() /= a_double;
  for (auto& element : this->getData()) {
    element /= a_double;
  }
  return (*this);
}

template <class GradientDataType>
PtWithGradient<GradientDataType>::PtWithGradient(void) {
  this->getPt() = Pt::fromScalarConstant(0.0);
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    this->getData()[d] = GradientDataType(0.0);
  }
}
template <class GradientDataType>
PtWithGradient<GradientDataType>::PtWithGradient(const Pt& a_pt)
    : PtWithDataCommon<PtWithGradient<GradientDataType>,
                       std::array<GradientDataType, 3>>(a_pt) {
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    this->getData()[d] = GradientDataType(0.0);
  }
}

template <class GradientDataType>
PtWithGradient<GradientDataType>& PtWithGradient<GradientDataType>::operator-(
    void) {
  this->getPt() = -this->getPt();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    this->getData()[d] = -this->getData()[d];
  }
  return (*this);
}

template <class GradientDataType>
PtWithGradient<GradientDataType>& PtWithGradient<GradientDataType>::operator=(
    const Pt& a_pt) {
  this->getPt() = a_pt;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    this->getData()[d] = GradientDataType(0.0);
  }
  return (*this);
}

template <class GradientDataType>
PtWithGradient<GradientDataType>& PtWithGradient<GradientDataType>::operator=(
    const PtWithGradient<GradientDataType>& a_pt) {
  this->getPt() = a_pt.getPt();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    this->getData()[d] = a_pt.getData()[d];
  }
  return (*this);
}

template <class GradientDataType>
PtWithGradient<GradientDataType>& PtWithGradient<GradientDataType>::operator+=(
    const PtWithGradient<GradientDataType>& a_other_pt) {
  this->getPt() += a_other_pt.getPt();
  const auto& other_data = a_other_pt.getData();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    this->getData()[d] += other_data[d];
  }
  return (*this);
}

template <class GradientDataType>
const PtWithGradient<GradientDataType> operator*(
    const double a_rhs, const PtWithGradient<GradientDataType>& a_pt) {
  auto pt_with_grad = PtWithGradient<GradientDataType>();
  pt_with_grad.getPt() = a_rhs * a_pt.getPt();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_with_grad.getData()[d] = a_rhs * a_pt.getData()[d];
  }
  return pt_with_grad;
}

template <class GradientDataType>
const PtWithGradient<GradientDataType> operator*(
    const PtWithGradient<GradientDataType>& a_pt, const double a_rhs) {
  auto pt_with_grad = PtWithGradient<GradientDataType>();
  pt_with_grad.getPt() = a_rhs * a_pt.getPt();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_with_grad.getData()[d] = a_rhs * a_pt.getData()[d];
  }
  return pt_with_grad;
}

template <class GradientDataType>
const PtWithGradient<GradientDataType> operator/(
    const PtWithGradient<GradientDataType>& a_pt, const double a_rhs) {
  auto pt_with_grad = PtWithGradient<GradientDataType>();
  pt_with_grad.getPt() = a_pt.getPt() / a_rhs;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_with_grad.getData()[d] = a_pt.getData()[d] / a_rhs;
  }
  return pt_with_grad;
}

template <class GradientDataType>
const PtWithGradient<GradientDataType> operator+(
    const PtWithGradient<GradientDataType>& a_pt1,
    const PtWithGradient<GradientDataType>& a_pt2) {
  auto pt_with_grad = PtWithGradient<GradientDataType>();
  pt_with_grad.getPt() = a_pt1.getPt() + a_pt2.getPt();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_with_grad.getData()[d] = a_pt1.getData()[d] + a_pt2.getData()[d];
  }
  return pt_with_grad;
}

template <class GradientDataType>
const PtWithGradient<GradientDataType> operator-(
    const PtWithGradient<GradientDataType>& a_pt1,
    const PtWithGradient<GradientDataType>& a_pt2) {
  auto pt_with_grad = PtWithGradient<GradientDataType>();
  pt_with_grad.getPt() = a_pt1.getPt() - a_pt2.getPt();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_with_grad.getData()[d] = a_pt1.getData()[d] - a_pt2.getData()[d];
  }
  return pt_with_grad;
}

template <class GradientDataType>
const PtWithGradient<GradientDataType> operator-(
    const PtWithGradient<GradientDataType>& a_pt) {
  auto pt_with_grad = PtWithGradient<GradientDataType>();
  pt_with_grad.getPt() = -a_pt.getPt();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_with_grad.getData()[d] = -a_pt.getData()[d];
  }
  return pt_with_grad;
}

template <class FunctorType, UnsignedIndex_t kArrayLength>
std::ostream& operator<<(
    std::ostream& out,
    const PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>&
        a_pt_with_data) {
  out << "Location: " << a_pt_with_data.getPt();
  out << "   Data: ";
  const auto& data = a_pt_with_data.getData();
  for (UnsignedIndex_t n = 0; n < data.size(); ++n) {
    out << data[n] << " ";
  }
  out << std::endl;
  return out;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_PT_WITH_DATA_TPP_
