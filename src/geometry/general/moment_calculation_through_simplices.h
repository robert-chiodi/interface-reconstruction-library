// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_H_
#define SRC_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_H_

#include "src/geometry/general/pt.h"
#include "src/moments/volume.h"
#include "src/moments/volume_moments.h"
#include "src/moments/volume_moments_and_doubles.h"
#include "src/moments/volume_moments_and_normal.h"

namespace IRL {

template <class GeometryType, class CalculationFunctor>
auto calculateMoments(const GeometryType& a_geometry,
                      CalculationFunctor a_moment_accumulator) ->
    typename CalculationFunctor::ReturnType;

class Volume3D_Functor {
 public:
  using ReturnType = Volume;

  Volume3D_Functor(void) : volume_m(Volume::fromScalarConstant(0.0)) {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  Volume volume_m;
};

class VolumeMoments3D_Functor {
 public:
  using ReturnType = VolumeMoments;

  VolumeMoments3D_Functor(void)
      : volume_moments_m(VolumeMoments::fromScalarConstant(0.0)) {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  VolumeMoments volume_moments_m;
};

class Centroid3D_Functor {
 public:
  using ReturnType = Pt;

  Centroid3D_Functor(void) : volume_moments_functor_m() {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  VolumeMoments3D_Functor volume_moments_functor_m;
};

template <UnsignedIndex_t kArrayLength>
class VolumeMomentsAndDoubles3D_Functor {
 public:
  using ReturnType = VolumeMomentsAndDoubles<kArrayLength>;

  VolumeMomentsAndDoubles3D_Functor(void)
      : volume_moments_and_doubles_m(
            VolumeMomentsAndDoubles<kArrayLength>::fromScalarConstant(0.0,
                                                                      0.0)) {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  ReturnType volume_moments_and_doubles_m;
};

class Volume2D_Functor {
 public:
  using ReturnType = Volume;

  Volume2D_Functor(void) : volume_m(Volume::fromScalarConstant(0.0)) {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  Volume volume_m;
};

class VolumeMoments2D_Functor {
 public:
  using ReturnType = VolumeMoments;

  VolumeMoments2D_Functor(void)
      : volume_moments_m(VolumeMoments::fromScalarConstant(0.0)) {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  VolumeMoments volume_moments_m;
};

class Centroid2D_Functor {
 public:
  using ReturnType = Pt;

  Centroid2D_Functor(void) : volume_moments_functor_m() {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  VolumeMoments2D_Functor volume_moments_functor_m;
};

class VolumeMomentsAndNormal2D_Functor {
 public:
  using ReturnType = VolumeMomentsAndNormal;

  VolumeMomentsAndNormal2D_Functor(void)
      : volume_moments_functor_m(), normal_m(Normal::fromScalarConstant(0.0)) {}

  template <class SimplexType>
  inline void operator()(const SimplexType& a_simplex);

  inline ReturnType getMoments(void) const;

 private:
  VolumeMoments2D_Functor volume_moments_functor_m;
  Normal normal_m;
};

}  // namespace IRL

#include "src/geometry/general/moment_calculation_through_simplices.tpp"

#endif  // SRC_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_H_
