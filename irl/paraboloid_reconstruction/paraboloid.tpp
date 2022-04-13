#ifndef IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_
#define IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_

namespace IRL {

inline Paraboloid::Paraboloid(const Pt& a_datum,
                              const ReferenceFrame& a_reference_frame,
                              const double a_coef_a, const double a_coef_b)
    : datum_m(a_datum),
      frame_m(a_reference_frame),
      paraboloid_m(std::array<double, 2>({a_coef_a, a_coef_b})) {
  assert(frame_m.isOrthonormalBasis());
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

inline Normal getParaboloidSurfaceNormal(const AlignedParaboloid& a_paraboloid,
                                         const Pt& a_pt) {
  return {2.0 * a_paraboloid.a() * a_pt[0], 2.0 * a_paraboloid.b() * a_pt[1],
          1.0};
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
  if (std::fabs(c) < 5.0 * DBL_EPSILON) {
    return a_starting_pt;
  } else {
    const auto solutions = solveQuadratic(a, b, c);
    if (solutions.size() == 0) {
      std::cout << a_line << a_starting_pt << std::endl;
    }
    assert(solutions.size() > 0);
    if (solutions.size() == 1) {
      assert(solutions[0] >= 0.0);
      return a_starting_pt + a_line * solutions[0];
    } else {
      assert(std::max(solutions[0], solutions[1]) >= 0.0);
      const double distance_along_line =
          solutions[0] > 0.0 && solutions[1] > 0.0
              ? std::min(solutions[0], solutions[1])
              : std::max(solutions[0], solutions[1]);
      return a_starting_pt + a_line * distance_along_line;
    }
  }
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_TPP_
