/*
 * geometry_classes.h
 *
 *  Created on: Jul 4, 2018
 *      Author: Robert Chiodi
 */

#ifndef SRC_GEOMETRY_CLASSES_H_
#define SRC_GEOMETRY_CLASSES_H_

#include <cmath>

namespace R2P {

// First class definitions will appear.
// Below this, all inlined function definitions will be given.

// A simple 3D point
class Pt {
 public:
  double x_m = {0.0}, y_m = {0.0}, z_m = {0.0};

  Pt(void) = default;
  Pt(const double a_x, const double a_y, const double a_z)
      : x_m(a_x), y_m(a_y), z_m(a_z) {}
  explicit Pt(const double* a_dim)
      : x_m(a_dim[0]), y_m(a_dim[1]), z_m(a_dim[2]) {}

  void setLoc(const double a_x, const double a_y, const double a_z) {
    x_m = a_x;
    y_m = a_y;
    z_m = a_z;
  }

  Pt& operator+=(const Pt& a_rhs) {
    (*this).x_m += a_rhs.x_m;
    (*this).y_m += a_rhs.y_m;
    (*this).z_m += a_rhs.z_m;
    return (*this);
  }

  Pt& operator/(const double a_rhs) {
    (*this).x_m /= a_rhs;
    (*this).y_m /= a_rhs;
    (*this).z_m /= a_rhs;
    return (*this);
  }

  ~Pt(void) = default;
};

// Overload plus and minus operators
inline Pt operator+(const Pt& a_pt1, const Pt& a_pt2);
inline Pt operator-(const Pt& a_pt1, const Pt& a_pt2);
inline Pt operator*(const double& a_multiplier, const Pt& a_pt2);
inline Pt operator*(const Pt& a_pt, const double& a_multiplier);

// A rectangular cuboid, such as a rectilinear cell
class RectangularCuboid {
 public:
  Pt vertex_m[8];

  // Constructor
  RectangularCuboid(void)
      : vertex_m{Pt(0.5, -0.5, -0.5), Pt(0.5, -0.5, 0.5),   Pt(0.5, 0.5, 0.5),
                 Pt(0.5, 0.5, -0.5),  Pt(-0.5, -0.5, -0.5), Pt(-0.5, -0.5, 0.5),
                 Pt(-0.5, 0.5, 0.5),  Pt(-0.5, 0.5, -0.5)} {}
  RectangularCuboid(const Pt& a_min_point, const Pt& a_max_point)
      : vertex_m{
            Pt(a_max_point.x_m, a_min_point.y_m, a_min_point.z_m),
            Pt(a_max_point.x_m, a_min_point.y_m, a_max_point.z_m),
            Pt(a_max_point.x_m, a_max_point.y_m, a_max_point.z_m),
            Pt(a_max_point.x_m, a_max_point.y_m, a_min_point.z_m),
            Pt(a_min_point.x_m, a_min_point.y_m, a_min_point.z_m),
            Pt(a_min_point.x_m, a_min_point.y_m, a_max_point.z_m),
            Pt(a_min_point.x_m, a_max_point.y_m, a_max_point.z_m),
            Pt(a_min_point.x_m, a_max_point.y_m, a_min_point.z_m),
        } {}

  // Shift in x/y/z
  inline void shiftInX(const double a_shift_dist);
  inline void shiftInY(const double a_shift_dist);
  inline void shiftInZ(const double a_shift_dist);
  inline void shift(const double a_x_shift, const double a_y_shift,
                    const double a_z_shift);

  // Calculate and return volume
  inline double volume(void) const;
  // Calculate and return centroid
  inline Pt centroid(void) const;

  // Default destructor
  ~RectangularCuboid(void) = default;
};

// Unit cube
class UnitCube : public RectangularCuboid {
 public:
  const Pt vertex_m[8] = {Pt(0.5, -0.5, -0.5),  Pt(0.5, -0.5, 0.5),
                          Pt(0.5, 0.5, 0.5),    Pt(0.5, 0.5, -0.5),
                          Pt(-0.5, -0.5, -0.5), Pt(-0.5, -0.5, 0.5),
                          Pt(-0.5, 0.5, 0.5),   Pt(-0.5, 0.5, -0.5)};

  // Default constructor
  UnitCube(void) = default;

  // Default destructor
  ~UnitCube(void) = default;
};

// A plane n \cdot x - d = 0
// Normal points liquid to gas
class Plane {
 public:
  double normal_m[3] = {0.0};
  double distance_m = {0.0};

  // Default constructor
  Plane(void) = default;

  // Default destructor
  ~Plane(void) = default;
};

// Storage for zeroeth and
// first order moments of liquid and gas
class GeometricMoments {
 public:
  double volume_m = {0.0};
  Pt centroid_m;

  // Default constructor and initialized constructor
  GeometricMoments(void) = default;
  GeometricMoments(const double a_volume, const Pt a_point)
      : volume_m(a_volume), centroid_m(a_point) {}

  GeometricMoments& operator+=(const GeometricMoments& a_rhs) {
    (*this).volume_m += a_rhs.volume_m;
    (*this).centroid_m += a_rhs.centroid_m;
    return (*this);
  }

  // Default destructor
  ~GeometricMoments(void) = default;
};

// Storage for zeroeth and
// first order moments of liquid and gas
class PhaseMoments {
 public:
  GeometricMoments liquid_m, gas_m;

  // Default constructor and initialized constructor
  PhaseMoments(void) = default;
  PhaseMoments(const GeometricMoments a_liquid_moments,
               const GeometricMoments a_gas_moments)
      : liquid_m(a_liquid_moments), gas_m(a_gas_moments) {}

  PhaseMoments& operator+=(const PhaseMoments& a_rhs) {
    (*this).liquid_m += a_rhs.liquid_m;
    (*this).gas_m += a_rhs.gas_m;
    return (*this);
  }

  // Default destructor
  ~PhaseMoments(void) = default;
};

// A multi-planar reconstruction for a cell
class Reconstruction {
 public:
  int number_of_interfaces_m = {0};
  Plane planes_m[2];
  double flip_cut_m = {1.0};

  // Default constructor
  Reconstruction(void) = default;

  // Boolean for if interface cutting needs to be flipped
  bool isFlipped(void) const { return flip_cut_m < 0.0; }

  // Default destructor
  ~Reconstruction(void) = default;
};

// A convex Tetrahedron
class Tet {
 public:
  Pt vertex_m[4];

  // Default constructor and initialized
  Tet(void) = default;
  Tet(const Pt& a_pt0, const Pt& a_pt1, const Pt& a_pt2, const Pt& a_pt3)
      : vertex_m{a_pt0, a_pt1, a_pt2, a_pt3} {}

  // Calculte volume and return
  inline double volume(void) const;

  // Calculate centroid and teturn
  inline Pt centroid(void) const;

  // Calculate volume and centroid
  inline GeometricMoments moments() const;

  // Default destructor
  ~Tet(void) = default;
};

//******************************************************************* //
//     Inlined function definitions placed below this.
//******************************************************************* //

inline Pt operator+(const Pt& a_pt1, const Pt& a_pt2) {
  return Pt(a_pt1.x_m + a_pt2.x_m, a_pt1.y_m + a_pt2.y_m,
            a_pt1.z_m + a_pt2.z_m);
}

inline Pt operator-(const Pt& a_pt1, const Pt& a_pt2) {
  return Pt(a_pt1.x_m - a_pt2.x_m, a_pt1.y_m - a_pt2.y_m,
            a_pt1.z_m - a_pt2.z_m);
}

inline Pt operator*(const double& a_multiplier, const Pt& a_pt) {
  return Pt(a_multiplier * a_pt.x_m, a_multiplier * a_pt.y_m,
            a_multiplier * a_pt.z_m);
}

inline Pt operator*(const Pt& a_pt, const double& a_multiplier) {
  return Pt(a_multiplier * a_pt.x_m, a_multiplier * a_pt.y_m,
            a_multiplier * a_pt.z_m);
}

inline void RectangularCuboid::shiftInX(const double a_shift_dist) {
  for (int v = 0; v < 8; ++v) {
    vertex_m[v].x_m = vertex_m[v].x_m + a_shift_dist;
  }
}

inline void RectangularCuboid::shiftInY(const double a_shift_dist) {
  for (int v = 0; v < 8; ++v) {
    vertex_m[v].y_m = vertex_m[v].y_m + a_shift_dist;
  }
}

inline void RectangularCuboid::shiftInZ(const double a_shift_dist) {
  for (int v = 0; v < 8; ++v) {
    vertex_m[v].z_m = vertex_m[v].z_m + a_shift_dist;
  }
}

inline void RectangularCuboid::shift(const double a_x_shift,
                                     const double a_y_shift,
                                     const double a_z_shift) {
  (*this).shiftInX(a_x_shift);
  (*this).shiftInY(a_y_shift);
  (*this).shiftInZ(a_z_shift);
}

// Calculate volume and return
inline double RectangularCuboid::volume(void) const {
  return (vertex_m[2].x_m - vertex_m[4].x_m) *
         (vertex_m[2].y_m - vertex_m[4].y_m) *
         (vertex_m[2].z_m - vertex_m[4].z_m);
}

inline Pt RectangularCuboid::centroid(void) const {
  return Pt(0.5 * (vertex_m[2].x_m + vertex_m[4].x_m),
            0.5 * (vertex_m[2].y_m + vertex_m[4].y_m),
            0.5 * (vertex_m[2].z_m + vertex_m[4].z_m));
}

// Calculte volume and return
inline double Tet::volume(void) const {
  Pt a = vertex_m[0] - vertex_m[3];
  Pt b = vertex_m[1] - vertex_m[3];
  Pt c = vertex_m[2] - vertex_m[3];

  return std::fabs(a.x_m * (b.y_m * c.z_m - c.y_m * b.z_m) -
                   a.y_m * (b.x_m * c.z_m - c.x_m * b.z_m) +
                   a.z_m * (b.x_m * c.y_m - c.x_m * b.y_m)) /
         6.0;
}

inline Pt Tet::centroid(void) const {
  return 0.25 * (vertex_m[0] + vertex_m[1] + vertex_m[2] + vertex_m[3]);
}

// Calculate volume and centroid, return through reference
inline GeometricMoments Tet::moments() const {
  return GeometricMoments((*this).volume(), (*this).centroid());
}

}  // namespace R2P

#endif /* SRC_GEOMETRY_CLASSES_H_ */
