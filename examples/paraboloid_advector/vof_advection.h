// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_

#include <string>

#include "irl/geometry/general/pt.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

#include "examples/paraboloid_advector/data.h"

void resetCentroids(const Data<IRL::LocalizedParaboloidLink<double>>&
                        a_link_localized_paraboloid,
                    Data<IRL::Pt>* a_liquid_centroid,
                    Data<IRL::Pt>* a_gas_centroid);

std::array<int, 3> getIndexFromTag(const BasicMesh& a_mesh,
                                   const IRL::UnsignedIndex_t a_tag);

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid);

void advectVOF(
    const std::string& a_advection_method,
    const std::string& a_reconstruction_method, const double a_dt,
    const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
    Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
    Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
    Data<IRL::Pt>* a_gas_centroid, Data<IRL::Paraboloid>* a_interface);

struct Split {
  static void advectVOF(
      const std::string& a_reconstruction_method, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
      Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
      Data<IRL::Pt>* a_gas_centroid, Data<IRL::Paraboloid>* a_interface);
};

struct FullLagrangian {
  static void advectVOF(
      const std::string& a_reconstruction_method, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
      Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
      Data<IRL::Pt>* a_gas_centroid);
};

struct SemiLagrangian {
  static void advectVOF(
      const std::string& a_reconstruction_method, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
      Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
      Data<IRL::Pt>* a_gas_centroid);
};

struct SemiLagrangianCorrected {
  static void advectVOF(
      const std::string& a_reconstruction_method, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::LocalizedParaboloidLink<double>>* a_link_localized_paraboloid,
      Data<double>* a_liquid_volume_fraction, Data<IRL::Pt>* a_liquid_centroid,
      Data<IRL::Pt>* a_gas_centroid);
};

inline IRL::Vec3<double> getVelocity(const IRL::Pt& a_location,
                                     const Data<double>& a_U,
                                     const Data<double>& a_V,
                                     const Data<double>& a_W);

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

inline IRL::Vec3<double> getVelocity(const IRL::Pt& a_location,
                                     const Data<double>& a_U,
                                     const Data<double>& a_V,
                                     const Data<double>& a_W) {
  return IRL::Vec3<double>(a_U.interpolate(a_location),
                           a_V.interpolate(a_location),
                           a_W.interpolate(a_location));
}

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_