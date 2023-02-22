#ifndef IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_
#define IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_

namespace IRL {

template <class ScalarType>
inline ParaboloidBase<ScalarType>::ParaboloidBase(void)
    : datum_m(),
      frame_m(),
      paraboloid_m(),
      place_infinite_shortcut_m({false, false}) {}

template <class ScalarType>
inline ParaboloidBase<ScalarType>::ParaboloidBase(
    const PtBase<ScalarType>& a_datum,
    const ReferenceFrameBase<ScalarType>& a_reference_frame,
    const ScalarType a_coef_a, const ScalarType a_coef_b)
    : datum_m(a_datum),
      frame_m(a_reference_frame),
      paraboloid_m(std::array<ScalarType, 2>({a_coef_a, a_coef_b})),
      place_infinite_shortcut_m({false, false}) {
  // assert(frame_m.isOrthonormalBasis());
}

template <class ScalarType>
inline ParaboloidBase<ScalarType> ParaboloidBase<ScalarType>::createAlwaysAbove(
    void) {
  ParaboloidBase<ScalarType> par;
  par.markAsAlwaysAbove();
  return par;
}

template <class ScalarType>
inline ParaboloidBase<ScalarType> ParaboloidBase<ScalarType>::createAlwaysBelow(
    void) {
  ParaboloidBase<ScalarType> par;
  par.markAsAlwaysBelow();
  return par;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::setDatum(
    const PtBase<ScalarType>& a_datum) {
  datum_m = a_datum;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::setReferenceFrame(
    const ReferenceFrameBase<ScalarType>& a_reference_frame) {
  assert(a_reference_frame.isOrthonormalBasis());
  frame_m = a_reference_frame;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::setAlignedParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid) {
  paraboloid_m = a_aligned_paraboloid;
}

template <class ScalarType>
inline const PtBase<ScalarType>& ParaboloidBase<ScalarType>::getDatum(
    void) const {
  return datum_m;
}

template <class ScalarType>
inline const ReferenceFrameBase<ScalarType>&
ParaboloidBase<ScalarType>::getReferenceFrame(void) const {
  return frame_m;
}

template <class ScalarType>
inline const AlignedParaboloidBase<ScalarType>&
ParaboloidBase<ScalarType>::getAlignedParaboloid(void) const {
  return paraboloid_m;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::markAsRealReconstruction(void) {
  place_infinite_shortcut_m[0] = false;
  place_infinite_shortcut_m[1] = false;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::markAsAlwaysAbove(void) {
  place_infinite_shortcut_m[0] = true;
  place_infinite_shortcut_m[1] = false;
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::markAsAlwaysBelow(void) {
  place_infinite_shortcut_m[0] = false;
  place_infinite_shortcut_m[1] = true;
}

template <class ScalarType>
inline bool ParaboloidBase<ScalarType>::isAlwaysAbove(void) const {
  return place_infinite_shortcut_m[0];
}

template <class ScalarType>
inline bool ParaboloidBase<ScalarType>::isAlwaysBelow(void) const {
  return place_infinite_shortcut_m[1];
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::serialize(ByteBuffer* a_buffer) const {
  datum_m.serialize(a_buffer);
  frame_m[0].serialize(a_buffer);
  frame_m[1].serialize(a_buffer);
  frame_m[2].serialize(a_buffer);
  paraboloid_m.serialize(a_buffer);
  const UnsignedIndex_t bool_to_int =
      (place_infinite_shortcut_m[0] ? 1 : 0) +
      2 * (place_infinite_shortcut_m[1] ? 1 : 0);
  a_buffer->pack(&bool_to_int, 1);
}

template <class ScalarType>
inline void ParaboloidBase<ScalarType>::unpackSerialized(ByteBuffer* a_buffer) {
  datum_m.unpackSerialized(a_buffer);
  frame_m[0].unpackSerialized(a_buffer);
  frame_m[1].unpackSerialized(a_buffer);
  frame_m[2].unpackSerialized(a_buffer);
  paraboloid_m.unpackSerialized(a_buffer);
  UnsignedIndex_t int_to_bool = 0;
  a_buffer->unpack(&int_to_bool, 1);
  place_infinite_shortcut_m[0] = int_to_bool % 2 == 1 ? true : false;
  place_infinite_shortcut_m[1] = int_to_bool / 2 == 1 ? true : false;
}
template <class ScalarType>
inline PtBase<ScalarType> conicCenter(
    const PlaneBase<ScalarType>& a_plane,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  const auto& face_normal = a_plane.normal();
  const auto& face_distance = a_plane.distance();
  return PtBase<ScalarType>(
      face_normal[0] / safelyTiny(static_cast<ScalarType>(2) *
                                  a_paraboloid.a() * face_normal[2]),
      face_normal[1] / safelyTiny(static_cast<ScalarType>(2) *
                                  a_paraboloid.b() * face_normal[2]),
      -face_normal[0] * face_normal[0] /
              safelyTiny(static_cast<ScalarType>(2) * a_paraboloid.a() *
                         face_normal[2] * face_normal[2]) -
          face_normal[1] * face_normal[1] /
              safelyTiny(static_cast<ScalarType>(2) * a_paraboloid.b() *
                         face_normal[2] * face_normal[2]) +
          face_distance / safelyTiny(face_normal[2]));
}

template <class ScalarType>
inline NormalBase<ScalarType> getParaboloidSurfaceNormal(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtBase<ScalarType>& a_pt) {
  return NormalBase<ScalarType>(
      static_cast<ScalarType>(2) * a_paraboloid.a() * a_pt[0],
      static_cast<ScalarType>(2) * a_paraboloid.b() * a_pt[1],
      static_cast<ScalarType>(1));
};

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient getParaboloidSurfaceNormalWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_pt) {
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(static_cast<ScalarType>(0)),
       B_grad = gradient_type(static_cast<ScalarType>(0));
  const auto& pt = a_pt.getPt();
  const auto& pt_grad = a_pt.getData();
  const auto surface_normal = PtBase<ScalarType>(
      static_cast<ScalarType>(2) * A * pt[0],
      static_cast<ScalarType>(2) * B * pt[1], static_cast<ScalarType>(1));
  auto surface_normal_withgrad = PtTypeWithGradient(surface_normal);
  auto& surface_normal_grad = surface_normal_withgrad.getData();
  surface_normal_grad[0] =
      static_cast<ScalarType>(2) * (A_grad * pt[0] + A * pt_grad[0]);
  surface_normal_grad[1] =
      static_cast<ScalarType>(2) * (B_grad * pt[1] + B * pt_grad[1]);
  // surface_normal_grad[2] = 0.0;
  return surface_normal_withgrad;
};

// Returns solution to quadratic equation solve.
// The smallest solution will always be first.
template <class ScalarType>
inline StackVector<ScalarType, 2> solveQuadraticBetween0And1(
    const ScalarType a, const ScalarType b, const ScalarType c) {
  ScalarType discriminant = b * b - static_cast<ScalarType>(4) * a * c;

  if (discriminant > static_cast<ScalarType>(0)) {
    if (a != static_cast<ScalarType>(0)) {
      /* First fast try in 32-bit precision */
      const ScalarType approx_discriminant = approxsqrt(discriminant);
      const ScalarType approx_q = -static_cast<ScalarType>(0.5) *
                                  (b + copysign(approx_discriminant, b));
      const ScalarType approx_sol1 = approx_q / safelyTiny(a);
      const ScalarType approx_sol2 = c / safelyTiny(approx_q);
      if ((approx_sol1 < -0.01 || approx_sol1 > 1.01) &&
          (approx_sol2 < -0.01 || approx_sol2 > 1.01)) {
        return StackVector<ScalarType, 2>();
      }

      /* Real calculation */
      discriminant = sqrt(discriminant);
      const ScalarType q =
          -static_cast<ScalarType>(0.5) * (b + copysign(discriminant, b));
      const ScalarType sol1 = q / safelyTiny(a);
      const ScalarType sol2 = c / safelyTiny(q);
      return sol1 < sol2 ? StackVector<ScalarType, 2>({sol1, sol2})
                         : StackVector<ScalarType, 2>({sol2, sol1});

    } else {
      return StackVector<ScalarType, 2>({-c / b});
    }
  }
  return StackVector<ScalarType, 2>();
};

// Returns solution to quadratic equation solve.
// The smallest solution will always be first.
template <class ScalarType>
inline StackVector<ScalarType, 2> solveQuadratic(const ScalarType a,
                                                 const ScalarType b,
                                                 const ScalarType c) {
  ScalarType discriminant = b * b - static_cast<ScalarType>(4) * a * c;
  // if (discriminant > static_cast<ScalarType>(0)) {
  //   discriminant = sqrt(discriminant);
  //   const ScalarType q =
  //       -(b + copysign(discriminant, b)) / static_cast<ScalarType>(2);
  //   const ScalarType sol1 = q / safelyTiny(a);
  //   const ScalarType sol2 = c / safelyTiny(q);
  //   // if (!isnan(sol1) && !isnan(sol2) &&
  //   //     !(fabs(q) < machine_epsilon<ScalarType>() &&
  //   //       fabs(a) < machine_epsilon<ScalarType>()) &&
  //   //     !(fabs(c) < machine_epsilon<ScalarType>() &&
  //   //       fabs(q) < machine_epsilon<ScalarType>())) {
  //   return sol1 < sol2 ? StackVector<ScalarType, 2>({sol1, sol2})
  //                      : StackVector<ScalarType, 2>({sol2, sol1});
  //   // }
  // }

  if (discriminant > static_cast<ScalarType>(0)) {
    if (a != static_cast<ScalarType>(0)) {
      discriminant = sqrt(discriminant);
      const ScalarType q =
          -static_cast<ScalarType>(0.5) * (b + copysign(discriminant, b));
      // if (b == static_cast<ScalarType>(0) && c ==
      // static_cast<ScalarType>(0))
      // {
      //   return StackVector<ScalarType, 2>(
      //       {static_cast<ScalarType>(0), static_cast<ScalarType>(0)});
      // } else if (q == static_cast<ScalarType>(0)) {
      //   return StackVector<ScalarType, 2>({static_cast<ScalarType>(0)});
      // }
      const ScalarType sol1 = q / safelyTiny(a);
      const ScalarType sol2 = c / safelyTiny(q);
      return sol1 < sol2 ? StackVector<ScalarType, 2>({sol1, sol2})
                         : StackVector<ScalarType, 2>({sol2, sol1});

    } else {
      return StackVector<ScalarType, 2>({-c / b});
    }
  }
  return StackVector<ScalarType, 2>();
};

// Returns solution to quadratic equation solve.
// The smallest solution will always be first.
template <class GradientType, class ScalarType>
inline StackVector<std::pair<ScalarType, GradientType>, 2>
solveQuadraticWithGradient(const ScalarType a, const ScalarType b,
                           const ScalarType c, const GradientType& a_grad,
                           const GradientType& b_grad,
                           const GradientType& c_grad) {
  using return_type = typename std::pair<ScalarType, GradientType>;
  const ScalarType discriminant = b * b - static_cast<ScalarType>(4) * a * c;
  if (discriminant > static_cast<ScalarType>(0)) {
    const auto d_grad = static_cast<ScalarType>(2) * b_grad * b -
                        static_cast<ScalarType>(4) * (a_grad * c + a * c_grad);
    if (a != static_cast<ScalarType>(0)) {
      const ScalarType sqrtd = sqrt(discriminant);
      const auto sqrtd_grad =
          d_grad / safelyEpsilon(sqrtd) / static_cast<ScalarType>(2);
      const ScalarType q =
          -(b + copysign(sqrtd, b)) / static_cast<ScalarType>(2);
      const auto q_grad =
          -(b_grad + sqrtd_grad * copysign(static_cast<ScalarType>(1), b)) /
          static_cast<ScalarType>(2);
      if (b == static_cast<ScalarType>(0) && c == static_cast<ScalarType>(0)) {
        const auto zero_return = std::pair<ScalarType, GradientType>(
            {static_cast<ScalarType>(0),
             GradientType(static_cast<ScalarType>(0))});
        return StackVector<return_type, 2>({zero_return, zero_return});
      } else if (q == static_cast<ScalarType>(0)) {
        const auto zero_return = std::pair<ScalarType, GradientType>(
            {static_cast<ScalarType>(0),
             GradientType(static_cast<ScalarType>(0))});
        return StackVector<return_type, 2>({zero_return});
      }
      const ScalarType sol1 = q / a;
      const ScalarType sol2 = c / q;
      const auto sol1_grad = (q_grad * a - q * a_grad) / safelyEpsilon(a * a);
      const auto sol2_grad = (c_grad * q - c * q_grad) / safelyEpsilon(q * q);
      const auto return_sol1 =
          std::pair<ScalarType, GradientType>({sol1, sol1_grad});
      const auto return_sol2 =
          std::pair<ScalarType, GradientType>({sol2, sol2_grad});
      return sol1 < sol2
                 ? StackVector<return_type, 2>({return_sol1, return_sol2})
                 : StackVector<return_type, 2>({return_sol2, return_sol1});

    } else {
      const ScalarType sol = -c / b;
      const auto sol_grad = (-c_grad * b + c * b_grad) / safelyEpsilon(b * b);
      const auto return_sol =
          std::pair<ScalarType, GradientType>({sol, sol_grad});
      return StackVector<return_type, 2>({return_sol});
    }
  }
  return StackVector<return_type, 2>();
};  // namespace IRL

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongLineOntoParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const ScalarType a = (a_paraboloid.a() * a_line[0] * a_line[0] +
                        a_paraboloid.b() * a_line[1] * a_line[1]);
  const ScalarType b = (a_line[2] +
                        static_cast<ScalarType>(2) * a_paraboloid.a() *
                            a_starting_pt[0] * a_line[0] +
                        static_cast<ScalarType>(2) * a_paraboloid.b() *
                            a_starting_pt[1] * a_line[1]);
  const ScalarType c = (a_starting_pt[2] +
                        a_paraboloid.a() * a_starting_pt[0] * a_starting_pt[0] +
                        a_paraboloid.b() * a_starting_pt[1] * a_starting_pt[1]);
  // check if starting point is on paraboloid (then solution = 0)
  // if (std::fabs(c) < 5.0 * DBL_EPSILON) {
  //   return a_starting_pt;
  // } else {
  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);
  if (solutions.size() == 0) {
    std::cout << "No solution found for projection on paraboloid" << a_line
              << a_starting_pt << std::endl;
  }
  assert(solutions.size() > 0);
  if (solutions.size() == 1) {
    return a_starting_pt + a_line * solutions[0];
  } else {
    if (abs(solutions[0]) < abs(solutions[1])) {
      return a_starting_pt + a_line * solutions[0];
    } else {
      return a_starting_pt + a_line * solutions[1];
    }
  }
  // }
}

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongHalfLineOntoParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const ScalarType a = (a_paraboloid.a() * a_line[0] * a_line[0] +
                        a_paraboloid.b() * a_line[1] * a_line[1]);
  const ScalarType b = (a_line[2] +
                        static_cast<ScalarType>(2) * a_paraboloid.a() *
                            a_starting_pt[0] * a_line[0] +
                        static_cast<ScalarType>(2) * a_paraboloid.b() *
                            a_starting_pt[1] * a_line[1]);
  const ScalarType c = (a_starting_pt[2] +
                        a_paraboloid.a() * a_starting_pt[0] * a_starting_pt[0] +
                        a_paraboloid.b() * a_starting_pt[1] * a_starting_pt[1]);
  // check if starting point is on paraboloid (then solution = 0)
  // if (std::fabs(c) < DBL_EPSILON) {
  //   return Pt(static_cast<ScalarType>(DBL_MAX),
  //             static_cast<ScalarType>(DBL_MAX),
  //             static_cast<ScalarType>(DBL_MAX));
  // } else {
  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);
  if (solutions.size() == 0) {
    return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                              static_cast<ScalarType>(DBL_MAX),
                              static_cast<ScalarType>(DBL_MAX));
  }
  // assert(solutions.size() > 0);
  if (solutions.size() == 1) {
    // assert(solutions[0] >= static_cast<ScalarType>(0));
    if (solutions[0] < machine_epsilon<ScalarType>()) {
      return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX));
    }
    return a_starting_pt + a_line * solutions[0];
  } else {
    // assert(maximum(solutions[0], solutions[1]) >=
    // static_cast<ScalarType>(0));
    if ((solutions[1] < static_cast<ScalarType>(0))) {
      return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX));

    } else {
      const ScalarType distance_along_line =
          solutions[0] > static_cast<ScalarType>(0)
              ? minimum(solutions[0], solutions[1])
              : maximum(solutions[0], solutions[1]);
      return a_starting_pt + a_line * distance_along_line;
    }
  }
  // }
}

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient projectPtAlongHalfLineOntoParaboloidWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_line, const PtTypeWithGradient& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  ScalarType EPSILON;
  if constexpr (std::is_same<ScalarType, Quad_t>::value) {
    EPSILON = FLT128_EPSILON;
  } else {
    EPSILON = DBL_EPSILON;
  }
  const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(static_cast<ScalarType>(0)),
       B_grad = gradient_type(static_cast<ScalarType>(0));
  A_grad.setGradA(static_cast<ScalarType>(1));
  B_grad.setGradB(static_cast<ScalarType>(1));
  const PtBase<ScalarType>& line = a_line.getPt();
  const auto& line_grad = a_line.getData();
  const PtBase<ScalarType>& starting_pt = a_starting_pt.getPt();
  const auto& starting_pt_grad = a_starting_pt.getData();
  const ScalarType a = (A * line[0] * line[0] + B * line[1] * line[1]);
  const ScalarType b =
      (line[2] + static_cast<ScalarType>(2) * A * starting_pt[0] * line[0] +
       static_cast<ScalarType>(2) * B * starting_pt[1] * line[1]);
  const ScalarType c = (starting_pt[2] + A * starting_pt[0] * starting_pt[0] +
                        B * starting_pt[1] * starting_pt[1]);
  // check if starting point is on paraboloid(then solution = 0)
  if (fabs(c) < static_cast<ScalarType>(5) * EPSILON) {
    return a_starting_pt;
  } else {
    const auto a_grad =
        A_grad * line[0] * line[0] + B_grad * line[1] * line[1] +
        static_cast<ScalarType>(2) * A * line_grad[0] * line[0] +
        static_cast<ScalarType>(2) * B * line_grad[1] * line[1];
    const auto b_grad =
        line_grad[2] +
        static_cast<ScalarType>(2) * A_grad * starting_pt[0] * line[0] +
        static_cast<ScalarType>(2) * A * starting_pt_grad[0] * line[0] +
        static_cast<ScalarType>(2) * A * starting_pt[0] * line_grad[0] +
        static_cast<ScalarType>(2) * B_grad * starting_pt[1] * line[1] +
        static_cast<ScalarType>(2) * B * starting_pt_grad[1] * line[1] +
        static_cast<ScalarType>(2) * B * starting_pt[1] * line_grad[1];
    const auto c_grad =
        starting_pt_grad[2] + A_grad * starting_pt[0] * starting_pt[0] +
        static_cast<ScalarType>(2) * A * starting_pt_grad[0] * starting_pt[0] +
        B_grad * starting_pt[1] * starting_pt[1] +
        static_cast<ScalarType>(2) * B * starting_pt_grad[1] * starting_pt[1];
    const auto solutions =
        solveQuadraticWithGradient(a, b, c, a_grad, b_grad, c_grad);
    if (solutions.size() == 0) {
      std::cout << "Solution not found on half-line: " << line << starting_pt
                << std::endl;
    }
    // assert(solutions.size() > 0);
    if (solutions.size() == 1) {
      const auto sol = solutions[0].first;
      const auto sol_grad = solutions[0].second;
      // assert(sol >= static_cast<ScalarType>(0));
      auto intersection =
          PtTypeWithGradient(PtBase<ScalarType>(starting_pt + sol * line));
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        intersection.getData()[d] =
            starting_pt_grad[d] + sol_grad * line[d] + sol * line_grad[d];
      }
      return intersection;
    } else {
      const auto sol0 = solutions[0].first;
      const auto sol1 = solutions[1].first;
      // assert(maximum(sol0, sol1) >= static_cast<ScalarType>(0));
      if (sol0 >= static_cast<ScalarType>(0)) {
        const auto sol0_grad = solutions[0].second;
        auto intersection =
            PtTypeWithGradient(PtBase<ScalarType>(starting_pt + sol0 * line));
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          intersection.getData()[d] =
              starting_pt_grad[d] + sol0_grad * line[d] + sol0 * line_grad[d];
        }
        return intersection;
      } else {
        const auto sol1_grad = solutions[1].second;
        auto intersection =
            PtTypeWithGradient(PtBase<ScalarType>(starting_pt + sol1 * line));
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          intersection.getData()[d] =
              starting_pt_grad[d] + sol1_grad * line[d] + sol1 * line_grad[d];
        }
        return intersection;
      }
    }
  }
}

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out, const ParaboloidBase<ScalarType>& a_paraboloid) {
  const auto& datum = a_paraboloid.getDatum();
  const auto& frame = a_paraboloid.getReferenceFrame();
  const auto& aligned_paraboloid = a_paraboloid.getAlignedParaboloid();

  out << "Datum: " << datum << '\n';
  out << "Frame: \n"
      << frame[0] << '\n'
      << frame[1] << '\n'
      << frame[2] << '\n';
  out << "Aligned Paraboloid: " << aligned_paraboloid << '\n';
  out << "is always above? " << a_paraboloid.isAlwaysAbove() << '\n';
  out << "is always below? " << a_paraboloid.isAlwaysBelow() << '\n';
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_
