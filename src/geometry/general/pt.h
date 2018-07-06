// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PT_H_
#define SRC_GEOMETRY_GENERAL_PT_H_

#include <array>
#include <cassert>
#include <iomanip>
#include <ostream>

#include "src/geometry/general/math_vector.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/expression_templates.h"
#include "src/helpers/mymath.h"
#include "src/helpers/serializer.h"
#include "src/parameters/defined_types.h"

namespace IRL {

class Pt;

inline std::ostream& operator<<(std::ostream& out, const Pt& a_pt);

/// \brief A point in 3D space.
class Pt : public Expr<Pt> {
 public:
  using value_type = double;
  /// \brief The default constructor, performs NO INITIALIZATION.
  Pt(void) = default;

  /// \brief A constructor that initializes the point given the x,y, and z
  /// locations.
  constexpr Pt(const double a_x, const double a_y, const double a_z);

  /// \brief Performs the interpolation on an edge between points with known
  /// distances to the plane to calculate the intersection point.
  ///
  /// This function takes two vertices of an edge and the corresponding signed
  /// distances to an intersecting plane. The point of intersection (where
  /// dist=0) is then calculated and returned.
  ///
  /// \param[in] a_point_1 First point on the edge.
  /// \param[in] a_dist_1 Corresponding signed distance from intersecting
  ///  plane for a_point_1.
  /// \param[in] a_point_2 Second point on the edge.
  /// \param[in] a_dist_2 Corresponding signed distance from intersecting
  ///  plane for a_point_2.
  static Pt fromEdgeIntersection(const Pt& a_pt_0, const double a_dist_0,
                                 const Pt& a_pt_1, const double a_dist_1);

  static Pt fromRawDoublePointer(const double* a_loc);

  static Pt fromArray(const std::array<double, 3>& a_loc);

  static constexpr Pt fromScalarConstant(const double a_value);

  static Pt fromVec3(const Vec3& a_vec);

  operator Vec3(void);

  /// \brief Assignment of point to a constant double.
  Pt& operator=(const double a_value);

  /// \brief Multiply point by constant value.
  Pt& operator*=(const double a_value);

  /// \brief Functions necessary to have consistent interface to PtWithData
  Pt& getPt(void);
  const Pt& getPt(void) const;

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  Pt(const Expr<E>& a_expr);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  Pt(Expr<E>&& a_expr);

  template <class E>
  Pt& operator=(const Expr<E>& a_expr);

  template <class E>
  Pt& operator=(Expr<E>&& a_expr);

  /// \brief Length of point vector.
  friend constexpr UnsignedIndex_t size(const Pt& x);

  /// \brief Return reference to `a_d` index of the point.
  double& operator[](const UnsignedIndex_t a_d);

  /// \brief Const version of `operator[]`.
  const double& operator[](const UnsignedIndex_t a_d) const;

  /// \brief Return value of `loc_m[0]` (x location).
  double& x(void);
  /// \brief Return value of `loc_m[1]` (y location).
  double& y(void);
  /// \brief Return value of `loc_m[2]` (z location).
  double& z(void);
  /// \brief Const version of `x()`.
  const double& x(void) const;
  /// \brief Const version of `y()`.
  const double& y(void) const;
  /// \brief Const version of `z()`.
  const double& z(void) const;

  /// \brief Return the index of the dimension with maximum magnitude.
  UnsignedIndex_t maxComponent(void) const;

  // \brief Unary minus operator
  Pt operator-(void);

  /// \brief Increment this point by another point.
  Pt& operator+=(const Pt& a_rhs);

  /// \brief Decrement this point by another point.
  Pt& operator-=(const Pt& a_rhs);

  /// \brief Divide each point location by the constant `a_rhs`.
  Pt& operator/(const double a_rhs);

  /// \brief Divide each point location by the constant `a_rhs`.
  Pt& operator/=(const double a_rhs);

  /// \brief Return size of the serialized point class in bytes.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize the point and store in the ByteBuffer
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the serialized point and store.
  void unpackSerialized(ByteBuffer* a_buffer);

  /// \brief Default destructor.
  ~Pt(void) = default;

 private:
  /// \brief A constructor that initializes the
  /// point given 3 elements in the array `a_loc`.
  explicit constexpr Pt(const double* a_loc);

  /// \brief A constructor that initializes the
  /// point given 3 elements in the array `a_loc`.
  explicit constexpr Pt(const std::array<double, 3>& a_loc);

  /// \brief Construct from a scalar constant value.
  explicit constexpr Pt(const double a_value);

  /// \brief Construct from a Vec3 MathVector.
  explicit Pt(const Vec3& a_vec);

  std::array<double, 3> loc_m;  ///< x,y,z (0,1,2) location of the point.
};

}  // namespace IRL

#include "src/geometry/general/pt.tpp"

#endif  // SRC_GEOMETRY_GENERAL_PT_H_
