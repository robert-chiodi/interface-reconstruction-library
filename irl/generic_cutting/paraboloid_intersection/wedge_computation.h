// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_WEDGE_COMPUTATION_H_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_WEDGE_COMPUTATION_H_

#include "irl/data_structures/small_vector.h"
#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/ellipse.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/parametrized_surface.h"

namespace IRL {

static constexpr double MINIMUM_EDGE_LENGTH = DBL_EPSILON;

template <class ReturnType>
ReturnType calculateTriangleCorrection(const AlignedParaboloid& a_paraboloid,
                                       const Pt& a_pt_0, const Pt& a_pt_1,
                                       const Pt& a_pt_2);

Normal getParaboloidSurfaceNormal(const AlignedParaboloid& a_paraboloid,
                                  const Pt& a_pt);

inline double signedDistance(const Pt& a_pt,
                             const AlignedParaboloid& a_paraboloid);

StackVector<double, 2> solveQuadratic(const double a, const double b,
                                      const double c);

Normal computeTangentVectorAtPoint(const AlignedParaboloid& a_paraboloid,
                                   const Plane& a_plane, const Pt& a_pt);

Normal computeAndCorrectTangentVectorAtPt(const AlignedParaboloid& a_paraboloid,
                                          const Plane& a_plane,
                                          const Pt& a_origin_pt,
                                          const Pt& a_end_pt,
                                          const Normal& a_end_tangent,
                                          const Pt& a_intersection_pt);

template <class ReturnType, class SurfaceOutputType = NoSurfaceOutput>
ReturnType bezierIntegrate(const AlignedParaboloid& a_paraboloid,
                           const Plane& a_plane, const Pt& a_pt_ref,
                           const Pt& a_pt_0, const Pt& a_pt_1,
                           const Normal& a_tangent_0, const Normal& a_tangent_1,
                           SurfaceOutputType* a_surface = NULL);

template <class ReturnType, class SurfaceOutputType = NoSurfaceOutput>
inline ReturnType computeWedgeCorrection(const AlignedParaboloid& a_paraboloid,
                                         const Plane& a_plane, const Pt& a_pt_0,
                                         const Pt& a_pt_1,
                                         const Pt& a_previous_edge,
                                         const Pt& a_next_edge,
                                         SurfaceOutputType* a_surface = NULL);

template <class ReturnType>
ReturnType computeV3Contribution(const AlignedParaboloid& a_paraboloid,
                                 const Pt& a_pt_0, const Pt& a_pt_1,
                                 const Pt& a_cp, const double a_weight);

static double v3Series[41][3] = {
    {2.09523809523809528832e-01, 3.80952380952380931234e-01,
     -5.71428571428571410729e-02},
    {-1.26984126984126897281e-02, 2.53968253968253954156e-01,
     -1.26984126984126897281e-02},
    {-2.02020202020201967985e-02, -1.50072150072150078959e-01,
     3.05916305916305968082e-02},
    {3.72960372960372960049e-02, 6.39360639360639360085e-02,
     -2.27328227328227328030e-02},
    {-4.04040404040404005359e-02, -1.24320124320124320016e-02,
     1.06560106560106560014e-02},
    {3.56871886283650976979e-02, -1.17007175830705235225e-02,
     -1.96404902287255220608e-03},
    {-2.81279092424603222033e-02, 1.93985580982484993873e-02,
     -2.53258952949355449144e-03},
    {2.05862249205902465843e-02, -1.90026691574679169883e-02,
     4.08205485604866413762e-03},
    {-1.42950332747075244816e-02, 1.55257315036558540822e-02,
     -4.04150970849045942934e-03},
    {9.54264472907628141796e-03, -1.15117618953936087789e-02,
     3.34612217772804624097e-03},
    {-6.17719348100287651143e-03, 8.02232919610763178797e-03,
     -2.51160614062754345907e-03},
    {3.90101474890349820060e-03, -5.35157257972321016154e-03,
     1.76920699049727420810e-03},
    {-2.41399989775078170628e-03, 3.45489893623047647497e-03,
     -1.19128627544748053046e-03},
    {1.46860565623520115952e-03, -2.17399892419069132657e-03,
     7.75354403547797083919e-04},
    {-8.80620601999301463869e-04, 1.33998761306499628715e-03,
     -4.91374024102316972994e-04},
    {5.21508692652322930310e-04, -8.11929701135352872333e-04,
     3.04769590179816240831e-04},
    {-3.05509568942368672697e-04, 4.84935823718045519101e-04,
     -1.85693554017880819258e-04},
    {1.77276395815434737553e-04, -2.86091263576557155047e-04,
     1.11457242528130960613e-04},
    {-1.02003016912458784591e-04, 1.66992280852464172270e-04,
     -6.60467997392183550414e-05},
    {5.82514823061649380099e-05, -9.65689972832239939063e-05,
     3.87058390159377495678e-05},
    {-3.30420824175868270054e-05, 5.53862690352201691383e-05,
     -2.24639917748133214606e-05},
    {1.86285795481852382328e-05, -3.15342725894044787469e-05,
     1.29264281372589721935e-05},
    {-1.04445437246298650441e-05, 1.78365273922045843719e-05,
     -7.38180846738453570066e-06},
    {5.82652104120046942685e-06, -1.00291738167552344136e-05,
     4.18683046937667869923e-06},
    {-3.23537647052159395237e-06, 5.60903039441516042084e-06,
     -2.36015114894640314994e-06},
    {1.78895249909732800395e-06, -3.12166207896178048349e-06,
     1.32305364506587308319e-06},
    {-9.85314706294302336548e-07, 1.72957799169836069021e-06,
     -7.37930243147500172352e-07},
    {5.40730333792878414157e-07, -9.54352573707662287295e-07,
     4.09677505966842900732e-07},
    {-2.95753356068542047445e-07, 5.24603798556049428080e-07,
     -2.26476468134631995398e-07},
    {1.61258810015876044985e-07, -2.87363379726446506705e-07,
     1.24710685182335891934e-07},
    {-8.76704778002015941588e-08, 1.56898500301773665938e-07,
     -6.84246579715181159030e-08},
    {4.75337224703033490158e-08, -8.54066883855229413404e-08,
     3.74166200173413135769e-08},
    {-2.57065245781141249776e-08, 4.63595284039086283402e-08,
     -2.03967889080998694981e-08},
    {1.38690173306780284126e-08, -2.50980391851122611262e-08,
     1.10865459080972638555e-08},
    {-7.46569873014004157925e-09, 1.35539622702173452016e-08,
     -6.00967114184549464555e-09},
    {4.01027350117916946203e-09, -7.30270649187992777560e-09,
     3.24937617487547770152e-09},
    {-2.14984854345419424283e-09, 3.92601604195155335685e-09,
     -1.75271787826727456715e-09},
    {1.15032202162643386740e-09, -2.10632478621065486152e-09,
     9.43297274700691519062e-10},
    {-6.14401522453570845869e-10, 1.12785551663862354372e-09,
     -5.06601052225360602018e-10},
    {3.27602053250195924866e-10, -6.02812432885591007361e-10,
     2.71528967029476965440e-10},
    {-1.74397233757762826095e-10, 3.21627843814674501789e-10,
     -1.45260052080825577373e-10}};
}  // namespace IRL

#include "irl/generic_cutting/paraboloid_intersection/wedge_computation.tpp"

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_WEDGE_COMPUTATION_H_
