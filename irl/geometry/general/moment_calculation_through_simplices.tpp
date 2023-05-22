// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_TPP_
#define IRL_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_TPP_

#include "irl/geometry/general/normal.h"
#include "irl/helpers/mymath.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class GeometryType, class CalculationFunctor>
auto calculateMoments(const GeometryType& a_geometry,
                      CalculationFunctor a_moment_accumulator) ->
    typename CalculationFunctor::ReturnType {
  for (UnsignedIndex_t s = 0;
       s < a_geometry.getNumberOfSimplicesInDecomposition(); ++s) {
    a_moment_accumulator(a_geometry.getSimplexFromDecomposition(s));
  }
  return a_moment_accumulator.getMoments();
}

template <UnsignedIndex_t ORDER>
template <class SimplexType>
void GeneralMoments3D_Functor<ORDER>::operator()(const SimplexType& a_simplex) {
  const auto& datum = a_simplex[3];
  datum_m = datum.getPt();
  const double sixv = scalarTripleProduct(a_simplex[0].getPt() - datum.getPt(),
                                          a_simplex[1].getPt() - datum.getPt(),
                                          a_simplex[2].getPt() - datum.getPt());

  UnsignedIndex_t prev_layer = 0;
  UnsignedIndex_t curr_layer = 1;

  S_m[0][0][prev_layer] = 1.0;
  D_m[0][0][prev_layer] = 1.0;
  C_m[0][0][prev_layer] = 1.0;
  moments_m[0] += sixv;

  const auto& v0 = a_simplex[0].getPt();
  const auto& v1 = a_simplex[1].getPt();
  const auto& v2 = a_simplex[2].getPt();

  // build up successive polynomial orders
  for (int corder = 1, m = 1; corder <= ORDER; ++corder) {
    for (int i = corder; i >= 0; --i) {
      for (int j = corder - i; j >= 0; --j, ++m) {
        const auto k = corder - i - j;
        C_m[i][j][curr_layer] = 0.0;
        D_m[i][j][curr_layer] = 0.0;
        S_m[i][j][curr_layer] = 0.0;
        if (i > 0) {
          C_m[i][j][curr_layer] +=
              (v2[0] - datum[0]) * C_m[i - 1][j][prev_layer];
          D_m[i][j][curr_layer] +=
              (v1[0] - datum[0]) * D_m[i - 1][j][prev_layer];
          S_m[i][j][curr_layer] +=
              (v0[0] - datum[0]) * S_m[i - 1][j][prev_layer];
        }
        if (j > 0) {
          C_m[i][j][curr_layer] +=
              (v2[1] - datum[1]) * C_m[i][j - 1][prev_layer];
          D_m[i][j][curr_layer] +=
              (v1[1] - datum[1]) * D_m[i][j - 1][prev_layer];
          S_m[i][j][curr_layer] +=
              (v0[1] - datum[1]) * S_m[i][j - 1][prev_layer];
        }
        if (k > 0) {
          C_m[i][j][curr_layer] += (v2[2] - datum_m[2]) * C_m[i][j][prev_layer];
          D_m[i][j][curr_layer] += (v1[2] - datum_m[2]) * D_m[i][j][prev_layer];
          S_m[i][j][curr_layer] += (v0[2] - datum_m[2]) * S_m[i][j][prev_layer];
        }
        D_m[i][j][curr_layer] += C_m[i][j][curr_layer];
        S_m[i][j][curr_layer] += D_m[i][j][curr_layer];
        moments_m[m] += sixv * S_m[i][j][curr_layer];
      }
    }
    prev_layer = 1 - prev_layer;
    curr_layer = 1 - curr_layer;
  }
}

template <UnsignedIndex_t ORDER>
inline typename GeneralMoments3D_Functor<ORDER>::ReturnType
GeneralMoments3D_Functor<ORDER>::getMoments(void) const {
  // Compute final moments (still is "datum" referenced space)
  auto mom = mom_array();
  UnsignedIndex_t prev_layer = 0;
  UnsignedIndex_t curr_layer = 1;

  mom[0] = moments_m[0] / 6.0;
  C_m[0][0][prev_layer] = 1.0;
  for (int corder = 1, m = 1; corder <= ORDER; ++corder) {
    for (int i = corder; i >= 0; --i) {
      for (int j = corder - i; j >= 0; --j, ++m) {
        const auto k = corder - i - j;
        C_m[i][j][curr_layer] = 0.0;
        if (i > 0) C_m[i][j][curr_layer] += C_m[i - 1][j][prev_layer];
        if (j > 0) C_m[i][j][curr_layer] += C_m[i][j - 1][prev_layer];
        if (k > 0) C_m[i][j][curr_layer] += C_m[i][j][prev_layer];
        mom[m] = moments_m[m] / (C_m[i][j][curr_layer] * (corder + 1) *
                                 (corder + 2) * (corder + 3));
      }
    }
    curr_layer = 1 - curr_layer;
    prev_layer = 1 - prev_layer;
  }

  // Shift moments back to global coordinate system

  // calculate and save Pascal's triangle
  C_m[0][0][0] = 1.0;
  for (int corder = 1, m = 1; corder <= ORDER; ++corder) {
    for (int i = corder; i >= 0; --i, ++m) {
      const auto j = corder - i;
      C_m[i][corder][0] = 1.0;
      if (i > 0 && j > 0) {  // Could be if i > corder?
        C_m[i][corder][0] = C_m[i][corder - 1][0] + C_m[i - 1][corder - 1][0];
      }
    }
  }

  // shift moments back to the original position using
  // \int_\Omega x^i y^j z^k d\vec r =
  // \int_\omega (x+\xi)^i (y+\eta)^j (z+\zeta)^k d\vec r =
  // \sum_{a,b,c=0}^{i,j,k} \binom{i}{a} \binom{j}{b} \binom{k}{c}
  // \xi^{i-a} \eta^{j-b} \zeta^{k-c} \int_\omega x^a y^b z^c d\vec r
  auto mom_shifted = GeneralMoments3D<ORDER>();
  mom_shifted[0] = mom[0];
  for (int corder = 1, m = 1; corder <= ORDER; ++corder) {
    for (int i = corder; i >= 0; --i)
      for (int j = corder - i; j >= 0; --j, ++m) {
        const auto k = corder - i - j;
        for (int mcorder = 0, mm = 0; mcorder <= corder; ++mcorder) {
          for (int mi = mcorder; mi >= 0; --mi)
            for (int mj = mcorder - mi; mj >= 0; --mj, ++mm) {
              const auto mk = mcorder - mi - mj;
              if (mi <= i && mj <= j && mk <= k) {
                mom_shifted[m] +=
                    C_m[mi][i][0] * C_m[mj][j][0] * C_m[mk][k][0] *
                    std::pow(datum_m[0], static_cast<double>(i - mi)) *
                    std::pow(datum_m[1], static_cast<double>(j - mj)) *
                    std::pow(datum_m[2], static_cast<double>(k - mk)) * mom[mm];
              }
            }
        }
      }
  }

  return mom_shifted;
}

template <class SimplexType>
void Volume3D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum = a_simplex[3];
  volume_m += scalarTripleProduct(a_simplex[0].getPt() - datum.getPt(),
                                  a_simplex[1].getPt() - datum.getPt(),
                                  a_simplex[2].getPt() - datum.getPt());
}

Volume3D_Functor::ReturnType Volume3D_Functor::getMoments(void) const {
  return volume_m / 6.0;
}

template <class SimplexType>
void VolumeMoments3D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum = a_simplex[3];
  auto six_times_tet_volume =
      scalarTripleProduct(a_simplex[0].getPt() - datum.getPt(),
                          a_simplex[1].getPt() - datum.getPt(),
                          a_simplex[2].getPt() - datum.getPt());
  volume_moments_m.volume() += six_times_tet_volume;
  volume_moments_m.centroid() +=
      six_times_tet_volume * (a_simplex[0].getPt() + a_simplex[1].getPt() +
                              a_simplex[2].getPt() + datum.getPt());
}

VolumeMoments3D_Functor::ReturnType VolumeMoments3D_Functor::getMoments(
    void) const {
  auto moments_copy = volume_moments_m;
  moments_copy /= 6.0;
  moments_copy.centroid() *= 0.25;
  return moments_copy;
}

template <class SimplexType>
void Centroid3D_Functor::operator()(const SimplexType& a_simplex) {
  volume_moments_functor_m(a_simplex);
}

Centroid3D_Functor::ReturnType Centroid3D_Functor::getMoments(void) const {
  auto volume_moments = volume_moments_functor_m.getMoments();
  volume_moments.normalizeByVolume();
  return volume_moments.centroid();
}

template <UnsignedIndex_t kArrayLength>
template <class SimplexType>
void VolumeMomentsAndDoubles3D_Functor<kArrayLength>::operator()(
    const SimplexType& a_simplex) {
  const auto& datum = a_simplex[3];
  auto six_times_tet_volume =
      scalarTripleProduct(a_simplex[0].getPt() - datum.getPt(),
                          a_simplex[1].getPt() - datum.getPt(),
                          a_simplex[2].getPt() - datum.getPt());
  volume_moments_and_doubles_m.volume() += six_times_tet_volume;
  volume_moments_and_doubles_m.centroid() +=
      six_times_tet_volume * (a_simplex[0].getPt() + a_simplex[1].getPt() +
                              a_simplex[2].getPt() + datum.getPt());
  for (UnsignedIndex_t n = 0; n < volume_moments_and_doubles_m.data().size();
       ++n) {
    volume_moments_and_doubles_m.data()[n] +=
        six_times_tet_volume *
        (a_simplex[0].getData()[n] + a_simplex[1].getData()[n] +
         a_simplex[2].getData()[n] + datum.getData()[n]);
  }
}

template <UnsignedIndex_t kArrayLength>
typename VolumeMomentsAndDoubles3D_Functor<kArrayLength>::ReturnType
VolumeMomentsAndDoubles3D_Functor<kArrayLength>::getMoments(void) const {
  auto moments_copy = volume_moments_and_doubles_m;
  moments_copy /= 6.0;
  moments_copy.centroid() *= 0.25;
  for (UnsignedIndex_t n = 0; n < moments_copy.data().size(); ++n) {
    moments_copy.data()[n] *= 0.25;
  }
  return moments_copy;
}

template <class SimplexType>
void Volume2D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum_pt = a_simplex[0];
  Normal edge_cross_product =
      Normal::fromPt(crossProduct(a_simplex[1].getPt() - datum_pt.getPt(),
                                  a_simplex[2].getPt() - datum_pt.getPt()));

  // Dot product imparts signed-ness to the polygon area.
  double twice_times_tri_volume = magnitude(edge_cross_product);
  edge_cross_product.normalize();
  volume_m +=
      twice_times_tri_volume *
      dotProduct(a_simplex.getPlaneOfExistence().normal(), edge_cross_product);
}

Volume2D_Functor::ReturnType Volume2D_Functor::getMoments(void) const {
  return 0.5 * volume_m;
}

template <class SimplexType>
void VolumeMoments2D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum_pt = a_simplex[0];
  Normal edge_cross_product =
      Normal::fromPt(crossProduct(a_simplex[1].getPt() - datum_pt.getPt(),
                                  a_simplex[2].getPt() - datum_pt.getPt()));

  // Dot product imparts signed-ness to the polygon area.
  double twice_times_tri_volume = magnitude(edge_cross_product);
  edge_cross_product.normalize();
  twice_times_tri_volume *=
      dotProduct(a_simplex.getPlaneOfExistence().normal(), edge_cross_product);

  volume_moments_m.volume() += twice_times_tri_volume;
  volume_moments_m.centroid() +=
      twice_times_tri_volume *
      (datum_pt.getPt() + a_simplex[1].getPt() + a_simplex[2].getPt());
}

VolumeMoments2D_Functor::ReturnType VolumeMoments2D_Functor::getMoments(
    void) const {
  auto moments_copy = volume_moments_m;
  moments_copy *= 0.5;
  moments_copy.centroid() /= 3.0;
  return moments_copy;
}

template <class SimplexType>
void Centroid2D_Functor::operator()(const SimplexType& a_simplex) {
  volume_moments_functor_m(a_simplex);
}

Centroid2D_Functor::ReturnType Centroid2D_Functor::getMoments(void) const {
  auto volume_moments = volume_moments_functor_m.getMoments();
  volume_moments.normalizeByVolume();
  return volume_moments.centroid();
}

template <class SimplexType>
void VolumeMomentsAndNormal2D_Functor::operator()(
    const SimplexType& a_simplex) {
  volume_moments_functor_m(a_simplex);
  normal_m = a_simplex.getPlaneOfExistence().normal();
}

VolumeMomentsAndNormal2D_Functor::ReturnType
VolumeMomentsAndNormal2D_Functor::getMoments(void) const {
  auto volume_moments = volume_moments_functor_m.getMoments();
  return {volume_moments, volume_moments.volume() * normal_m};
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_TPP_
