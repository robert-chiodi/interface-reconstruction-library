#ifndef IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_
#define IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_

namespace IRL {

inline Paraboloid::Paraboloid(void)
    : datum_m(),
      frame_m(),
      paraboloid_m(),
      place_infinite_shortcut_m({false, false}) {}

inline Paraboloid::Paraboloid(const Pt& a_datum,
                              const ReferenceFrame& a_reference_frame,
                              const double a_coef_a, const double a_coef_b)
    : datum_m(a_datum),
      frame_m(a_reference_frame),
      paraboloid_m(std::array<double, 2>({a_coef_a, a_coef_b})),
      place_infinite_shortcut_m({false, false}) {
  // assert(frame_m.isOrthonormalBasis());
}

inline Paraboloid Paraboloid::createAlwaysAbove(void) {
  Paraboloid par;
  par.markAsAlwaysAbove();
  return par;
}

inline Paraboloid Paraboloid::createAlwaysBelow(void) {
  Paraboloid par;
  par.markAsAlwaysBelow();
  return par;
}

inline void Paraboloid::setDatum(const Pt& a_datum) { datum_m = a_datum; }

inline void Paraboloid::setReferenceFrame(
    const ReferenceFrame& a_reference_frame) {
  assert(a_reference_frame.isOrthonormalBasis());
  frame_m = a_reference_frame;
}

inline void Paraboloid::setAlignedParaboloid(
    const AlignedParaboloid& a_aligned_paraboloid) {
  paraboloid_m = a_aligned_paraboloid;
}

inline const Pt& Paraboloid::getDatum(void) const { return datum_m; }

inline const ReferenceFrame& Paraboloid::getReferenceFrame(void) const {
  return frame_m;
}

inline const AlignedParaboloid& Paraboloid::getAlignedParaboloid(void) const {
  return paraboloid_m;
}

inline void Paraboloid::markAsRealReconstruction(void) {
  place_infinite_shortcut_m[0] = true;
  place_infinite_shortcut_m[1] = true;
}

inline void Paraboloid::markAsAlwaysAbove(void) {
  place_infinite_shortcut_m[0] = true;
  place_infinite_shortcut_m[1] = false;
}

inline void Paraboloid::markAsAlwaysBelow(void) {
  place_infinite_shortcut_m[0] = false;
  place_infinite_shortcut_m[1] = true;
}

inline bool Paraboloid::isAlwaysAbove(void) const {
  return place_infinite_shortcut_m[0];
}

inline bool Paraboloid::isAlwaysBelow(void) const {
  return place_infinite_shortcut_m[1];
}

inline Pt conicCenter(const Plane& a_plane,
                      const AlignedParaboloid& a_paraboloid) {
  const auto& face_normal = a_plane.normal();
  const auto& face_distance = a_plane.distance();
  return {face_normal[0] / safelyTiny(2.0 * a_paraboloid.a() * face_normal[2]),
          face_normal[1] / safelyTiny(2.0 * a_paraboloid.b() * face_normal[2]),
          -face_normal[0] * face_normal[0] /
                  safelyTiny(2.0 * a_paraboloid.a() * face_normal[2] *
                             face_normal[2]) -
              face_normal[1] * face_normal[1] /
                  safelyTiny(2.0 * a_paraboloid.b() * face_normal[2] *
                             face_normal[2]) +
              face_distance / safelyTiny(face_normal[2])};
}

inline Normal getParaboloidSurfaceNormal(const AlignedParaboloid& a_paraboloid,
                                         const Pt& a_pt) {
  return {2.0 * a_paraboloid.a() * a_pt[0], 2.0 * a_paraboloid.b() * a_pt[1],
          1.0};
};

template <class PtTypeWithGradient>
inline PtTypeWithGradient getParaboloidSurfaceNormalWithGradient(
    const AlignedParaboloid& a_paraboloid, const PtTypeWithGradient& a_pt) {
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  const double A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  const auto& pt = a_pt.getPt();
  const auto& pt_grad = a_pt.getData();
  const auto surface_normal = Pt(2.0 * A * pt[0], 2.0 * B * pt[1], 1.0);
  auto surface_normal_withgrad = PtTypeWithGradient(surface_normal);
  auto& surface_normal_grad = surface_normal_withgrad.getData();
  surface_normal_grad[0] = 2.0 * (A_grad * pt[0] + A * pt_grad[0]);
  surface_normal_grad[1] = 2.0 * (B_grad * pt[1] + B * pt_grad[1]);
  // surface_normal_grad[2] = 0.0;
  return surface_normal_withgrad;
};

// Returns solution to quadratic equation solve.
// The smallest solution will always be first.
inline StackVector<double, 2> solveQuadratic(const double a, const double b,
                                             const double c) {
  double discriminant = b * b - 4.0 * a * c;
  // By preventing discriminant = 0, we avoid cases with intersections tangent
  // to the paraboloid
  if (discriminant > 0.0) {
    if (a != 0.0) {
      discriminant = std::sqrt(discriminant);
      const double q = -0.5 * (b + std::copysign(discriminant, b));
      if (b == 0.0 && c == 0.0) {
        return StackVector<double, 2>({0.0, 0.0});
      } else if (q == 0.0) {
        return StackVector<double, 2>({0.0});
      }
      const double sol1 = q / a;
      const double sol2 = c / q;
      return sol1 < sol2 ? StackVector<double, 2>({sol1, sol2})
                         : StackVector<double, 2>({sol2, sol1});

    } else {
      return StackVector<double, 2>({-c / b});
    }
  }
  return StackVector<double, 2>();
};

// Returns solution to quadratic equation solve.
// The smallest solution will always be first.
template <class GradientType>
inline StackVector<std::pair<double, GradientType>, 2>
solveQuadraticWithGradient(const double a, const double b, const double c,
                           const GradientType& a_grad,
                           const GradientType& b_grad,
                           const GradientType& c_grad) {
  using return_type = typename std::pair<double, GradientType>;
  const double discriminant = b * b - 4.0 * a * c;
  if (discriminant > 0.0) {
    const auto d_grad = 2.0 * b_grad * b - 4.0 * (a_grad * c + a * c_grad);
    if (a != 0.0) {
      const double sqrtd = std::sqrt(discriminant);
      const auto sqrtd_grad = 0.5 * d_grad / safelyEpsilon(sqrtd);
      const double q = -0.5 * (b + std::copysign(sqrtd, b));
      const auto q_grad = -0.5 * (b_grad + sqrtd_grad * std::copysign(1.0, b));
      if (b == 0.0 && c == 0.0) {
        const auto zero_return =
            std::pair<double, GradientType>({0.0, GradientType(0.0)});
        return StackVector<return_type, 2>({zero_return, zero_return});
      } else if (q == 0.0) {
        const auto zero_return =
            std::pair<double, GradientType>({0.0, GradientType(0.0)});
        return StackVector<return_type, 2>({zero_return});
      }
      const double sol1 = q / a;
      const double sol2 = c / q;
      const auto sol1_grad = (q_grad * a - q * a_grad) / safelyEpsilon(a * a);
      const auto sol2_grad = (c_grad * q - c * q_grad) / safelyEpsilon(q * q);
      const auto return_sol1 =
          std::pair<double, GradientType>({sol1, sol1_grad});
      const auto return_sol2 =
          std::pair<double, GradientType>({sol2, sol2_grad});
      return sol1 < sol2
                 ? StackVector<return_type, 2>({return_sol1, return_sol2})
                 : StackVector<return_type, 2>({return_sol2, return_sol1});

    } else {
      const double sol = -c / b;
      const auto sol_grad = (-c_grad * b + c * b_grad) / safelyEpsilon(b * b);
      const auto return_sol = std::pair<double, GradientType>({sol, sol_grad});
      return StackVector<return_type, 2>({return_sol});
    }
  }
  return StackVector<return_type, 2>();
};  // namespace IRL

inline Pt projectPtAlongLineOntoParaboloid(
    const AlignedParaboloid& a_paraboloid, const Normal& a_line,
    const Pt& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const double a = (a_paraboloid.a() * a_line[0] * a_line[0] +
                    a_paraboloid.b() * a_line[1] * a_line[1]);
  const double b =
      (a_line[2] + 2.0 * a_paraboloid.a() * a_starting_pt[0] * a_line[0] +
       2.0 * a_paraboloid.b() * a_starting_pt[1] * a_line[1]);
  const double c = (a_starting_pt[2] +
                    a_paraboloid.a() * a_starting_pt[0] * a_starting_pt[0] +
                    a_paraboloid.b() * a_starting_pt[1] * a_starting_pt[1]);
  // check if starting point is on paraboloid (then solution = 0)
  // if (std::fabs(c) < 5.0 * DBL_EPSILON) {
  //   return a_starting_pt;
  // } else {
  const auto solutions = solveQuadratic(a, b, c);
  if (solutions.size() == 0) {
    std::cout << "No solution found for projection on paraboloid" << a_line
              << a_starting_pt << std::endl;
  }
  assert(solutions.size() > 0);
  if (solutions.size() == 1) {
    return a_starting_pt + a_line * solutions[0];
  } else {
    if (std::abs(solutions[0]) < std::abs(solutions[1])) {
      return a_starting_pt + a_line * solutions[0];
    } else {
      return a_starting_pt + a_line * solutions[1];
    }
  }
  // }
}

inline Pt projectPtAlongHalfLineOntoParaboloid(
    const AlignedParaboloid& a_paraboloid, const Normal& a_line,
    const Pt& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const double a = (a_paraboloid.a() * a_line[0] * a_line[0] +
                    a_paraboloid.b() * a_line[1] * a_line[1]);
  const double b =
      (a_line[2] + 2.0 * a_paraboloid.a() * a_starting_pt[0] * a_line[0] +
       2.0 * a_paraboloid.b() * a_starting_pt[1] * a_line[1]);
  const double c = (a_starting_pt[2] +
                    a_paraboloid.a() * a_starting_pt[0] * a_starting_pt[0] +
                    a_paraboloid.b() * a_starting_pt[1] * a_starting_pt[1]);
  // check if starting point is on paraboloid (then solution = 0)
  // if (std::fabs(c) < DBL_EPSILON) {
  //   return Pt(DBL_MAX, DBL_MAX, DBL_MAX);
  // } else {
  const auto solutions = solveQuadratic(a, b, c);
  if (solutions.size() == 0) {
    return Pt(DBL_MAX, DBL_MAX, DBL_MAX);
  }
  assert(solutions.size() > 0);
  if (solutions.size() == 1) {
    assert(solutions[0] >= 0.0);
    if (solutions[0] < 0.0) {
      return Pt(DBL_MAX, DBL_MAX, DBL_MAX);
    }
    return a_starting_pt + a_line * solutions[0];
  } else {
    assert(std::max(solutions[0], solutions[1]) >= 0.0);
    const double distance_along_line =
        solutions[0] > 0.0 && solutions[1] > 0.0
            ? std::min(solutions[0], solutions[1])
            : std::max(solutions[0], solutions[1]);
    if (distance_along_line < 0.0) {
      return Pt(DBL_MAX, DBL_MAX, DBL_MAX);
    }
    return a_starting_pt + a_line * distance_along_line;
  }
  // }
}

template <class PtTypeWithGradient>
inline PtTypeWithGradient projectPtAlongHalfLineOntoParaboloidWithGradient(
    const AlignedParaboloid& a_paraboloid, const PtTypeWithGradient& a_line,
    const PtTypeWithGradient& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  const double A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const Pt& line = a_line.getPt();
  const auto& line_grad = a_line.getData();
  const Pt& starting_pt = a_starting_pt.getPt();
  const auto& starting_pt_grad = a_starting_pt.getData();
  const double a = (A * line[0] * line[0] + B * line[1] * line[1]);
  const double b = (line[2] + 2.0 * A * starting_pt[0] * line[0] +
                    2.0 * B * starting_pt[1] * line[1]);
  const double c = (starting_pt[2] + A * starting_pt[0] * starting_pt[0] +
                    B * starting_pt[1] * starting_pt[1]);
  // check if starting point is on paraboloid(then solution = 0)
  if (std::fabs(c) < 5.0 * DBL_EPSILON) {
    return a_starting_pt;
  } else {
    const auto a_grad =
        A_grad * line[0] * line[0] + B_grad * line[1] * line[1] +
        2.0 * A * line_grad[0] * line[0] + 2.0 * B * line_grad[1] * line[1];
    const auto b_grad = line_grad[2] + 2.0 * A_grad * starting_pt[0] * line[0] +
                        2.0 * A * starting_pt_grad[0] * line[0] +
                        2.0 * A * starting_pt[0] * line_grad[0] +
                        2.0 * B_grad * starting_pt[1] * line[1] +
                        2.0 * B * starting_pt_grad[1] * line[1] +
                        2.0 * B * starting_pt[1] * line_grad[1];
    const auto c_grad = starting_pt_grad[2] +
                        A_grad * starting_pt[0] * starting_pt[0] +
                        2.0 * A * starting_pt_grad[0] * starting_pt[0] +
                        B_grad * starting_pt[1] * starting_pt[1] +
                        2.0 * B * starting_pt_grad[1] * starting_pt[1];
    const auto solutions =
        solveQuadraticWithGradient(a, b, c, a_grad, b_grad, c_grad);
    if (solutions.size() == 0) {
      std::cout << "Solution not found on half-line: " << line << starting_pt
                << std::endl;
    }
    assert(solutions.size() > 0);
    if (solutions.size() == 1) {
      const auto sol = solutions[0].first;
      const auto sol_grad = solutions[0].second;
      assert(sol >= 0.0);
      auto intersection = PtTypeWithGradient(Pt(starting_pt + sol * line));
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        intersection.getData()[d] =
            starting_pt_grad[d] + sol_grad * line[d] + sol * line_grad[d];
      }
      return intersection;
    } else {
      const auto sol0 = solutions[0].first;
      const auto sol1 = solutions[1].first;
      assert(std::max(sol0, sol1) >= 0.0);
      if (sol0 >= 0.0) {
        const auto sol0_grad = solutions[0].second;
        auto intersection = PtTypeWithGradient(Pt(starting_pt + sol0 * line));
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          intersection.getData()[d] =
              starting_pt_grad[d] + sol0_grad * line[d] + sol0 * line_grad[d];
        }
        return intersection;
      } else {
        const auto sol1_grad = solutions[1].second;
        auto intersection = PtTypeWithGradient(Pt(starting_pt + sol1 * line));
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          intersection.getData()[d] =
              starting_pt_grad[d] + sol1_grad * line[d] + sol1 * line_grad[d];
        }
        return intersection;
      }
    }
  }
}

inline std::ostream& operator<<(std::ostream& out,
                                const Paraboloid& a_paraboloid) {
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
