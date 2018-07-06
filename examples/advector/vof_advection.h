// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_ADVECTOR_VOF_ADVECTION_H_
#define EXAMPLES_ADVECTOR_VOF_ADVECTION_H_

#include <string>

#include "src/geometry/general/pt.h"
#include "src/planar_reconstruction/localized_separator_link.h"
#include "src/planar_reconstruction/localizer_link_from_localized_separator_link.h"

#include "examples/advector/data.h"

void resetCentroids(
    const Data<IRL::LocalizedSeparatorLink>& a_link_localized_separator,
    Data<IRL::Pt>* a_liquid_centroid, Data<IRL::Pt>* a_gas_centroid);

std::array<int, 3> getIndexFromTag(const BasicMesh& a_mesh,
                                   const IRL::UnsignedIndex_t a_tag);

void connectMesh(const BasicMesh& a_mesh,
                 Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator);

void advectVOF(const std::string& a_advection_method, const double a_dt,
               const Data<double>& a_U, const Data<double>& a_V,
               const Data<double>& a_W,
               Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
               Data<double>* a_liquid_volume_fraction,
               Data<IRL::Pt>* a_liquid_centroid, Data<IRL::Pt>* a_gas_centroid);

struct FullLagrangian {
  static void advectVOF(
      const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
      const Data<double>& a_W,
      Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
      Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
      Data<IRL::Pt>* a_gas_centroid);
};

struct SemiLagrangian {
  static void advectVOF(
      const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
      const Data<double>& a_W,
      Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
      Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
      Data<IRL::Pt>* a_gas_centroid);
};

struct SemiLagrangianCorrected {
  static void advectVOF(
      const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
      const Data<double>& a_W,
      Data<IRL::LocalizedSeparatorLink>* a_link_localized_separator,
      Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
      Data<IRL::Pt>* a_gas_centroid);
};

inline IRL::Vec3 getVelocity(const IRL::Pt& a_location, const Data<double>& a_U,
                             const Data<double>& a_V, const Data<double>& a_W);

inline IRL::Pt back_project_vertex(const IRL::Pt& a_initial_pt,
                                   const double a_dt, const Data<double>& a_U,
                                   const Data<double>& a_V,
                                   const Data<double>& a_W);

void correctCentroidLocation(Data<IRL::Pt>* a_liquid_centroid,
                             Data<IRL::Pt>* a_gas_centroid);

// ************************************************
//     Inlined functions below this
// ************************************************
// Spatial RK4 projection of a point in a velocity field..
inline IRL::Pt back_project_vertex(const IRL::Pt& a_initial_pt,
                                   const double a_dt, const Data<double>& a_U,
                                   const Data<double>& a_V,
                                   const Data<double>& a_W) {
  auto v1 = getVelocity(a_initial_pt, a_U, a_V, a_W);
  auto v2 = getVelocity(a_initial_pt + IRL::Pt::fromVec3(0.5 * a_dt * v1), a_U,
                        a_V, a_W);
  auto v3 = getVelocity(a_initial_pt + IRL::Pt::fromVec3(0.5 * a_dt * v2), a_U,
                        a_V, a_W);
  auto v4 =
      getVelocity(a_initial_pt + IRL::Pt::fromVec3(a_dt * v3), a_U, a_V, a_W);
  return a_initial_pt +
         IRL::Pt::fromVec3(a_dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0);
}

inline IRL::Vec3 getVelocity(const IRL::Pt& a_location, const Data<double>& a_U,
                             const Data<double>& a_V, const Data<double>& a_W) {
  return IRL::Vec3(a_U.interpolate(a_location), a_V.interpolate(a_location),
                   a_W.interpolate(a_location));
}

#endif  // EXAMPLES_ADVECTOR_VOF_ADVECTION_H_
