// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/paraboloid_reconstruction/c_paraboloid.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/helpers/mymath.h"
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"

#include <Eigen/Dense>
#include <iostream>

extern "C" {

void c_Paraboloid_new(c_Paraboloid* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::Paraboloid;
}

void c_Paraboloid_newFromObjectAllocationServer(
    c_Paraboloid* a_self, c_ObjServer_Paraboloid* a_object_allocation_server) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
}

void c_Paraboloid_delete(c_Paraboloid* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_Paraboloid_setDatum(c_Paraboloid* a_self, const double* a_datum) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setDatum(IRL::Pt::fromRawDoublePointer(a_datum));
}

void c_Paraboloid_setReferenceFrame(c_Paraboloid* a_self,
                                    const double* a_normal1,
                                    const double* a_normal2,
                                    const double* a_normal3) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setReferenceFrame(
      IRL::ReferenceFrame(IRL::Normal::fromRawDoublePointer(a_normal1),
                          IRL::Normal::fromRawDoublePointer(a_normal2),
                          IRL::Normal::fromRawDoublePointer(a_normal3)));
}

void c_Paraboloid_setAlignedParaboloid(c_Paraboloid* a_self,
                                       const double* a_coeff_a,
                                       const double* a_coeff_b) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setAlignedParaboloid(
      IRL::AlignedParaboloid({(*a_coeff_a), (*a_coeff_b)}));
}

void c_Paraboloid_setParaboloidFromPolygon(c_Paraboloid* a_self,
                                           c_Poly* a_polygon) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  IRL::Polygon polygon = (*a_polygon->obj_ptr);
  IRL::Normal normal_poly = polygon.calculateNormal();
  IRL::ReferenceFrame frame;
  IRL::UnsignedIndex_t largest_dir = 0;
  if (std::fabs(normal_poly[largest_dir]) < std::fabs(normal_poly[1]))
    largest_dir = 1;
  if (std::fabs(normal_poly[largest_dir]) < std::fabs(normal_poly[2]))
    largest_dir = 2;
  if (largest_dir == 0)
    frame[0] = crossProduct(normal_poly, IRL::Normal(0.0, 1.0, 0.0));
  else if (largest_dir == 1)
    frame[0] = crossProduct(normal_poly, IRL::Normal(0.0, 0.0, 1.0));
  else
    frame[0] = crossProduct(normal_poly, IRL::Normal(1.0, 0.0, 0.0));
  frame[0].normalize();
  frame[1] = crossProduct(normal_poly, frame[0]);
  frame[2] = normal_poly;
  IRL::Pt datum = polygon.calculateCentroid();

  a_self->obj_ptr->markAsRealReconstruction();
  a_self->obj_ptr->setDatum(datum);
  a_self->obj_ptr->setReferenceFrame(frame);
  a_self->obj_ptr->setAlignedParaboloid(
      IRL::AlignedParaboloid({1.0e-12, -1.0e-12}));
}

void c_Paraboloid_setParaboloidFromPlanarSep(c_Paraboloid* a_self,
                                             c_PlanarSep* a_plane,
                                             double* a_pt_ref) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  if (a_plane != nullptr && a_plane->obj_ptr != nullptr) {
    a_self->obj_ptr->markAsRealReconstruction();
    IRL::Plane plane = (*a_plane->obj_ptr)[0];
    IRL::Normal normal_plane = plane.normal();
    double distance_plane = plane.distance();
    if (IRL::squaredMagnitude(normal_plane) < DBL_MIN) {
      if (distance_plane < 0.0) {
        a_self->obj_ptr->markAsAlwaysBelow();
      } else {
        a_self->obj_ptr->markAsAlwaysAbove();
      }
    } else {
      IRL::Pt pt_ref = IRL::Pt::fromRawDoublePointer(a_pt_ref);
      IRL::ReferenceFrame frame;
      IRL::UnsignedIndex_t largest_dir = 0;
      if (std::fabs(normal_plane[largest_dir]) < std::fabs(normal_plane[1]))
        largest_dir = 1;
      if (std::fabs(normal_plane[largest_dir]) < std::fabs(normal_plane[2]))
        largest_dir = 2;
      if (largest_dir == 0)
        frame[0] = crossProduct(normal_plane, IRL::Normal(0.0, 1.0, 0.0));
      else if (largest_dir == 1)
        frame[0] = crossProduct(normal_plane, IRL::Normal(0.0, 0.0, 1.0));
      else
        frame[0] = crossProduct(normal_plane, IRL::Normal(1.0, 0.0, 0.0));
      frame[0].normalize();
      frame[1] = crossProduct(normal_plane, frame[0]);
      frame[2] = normal_plane;
      IRL::Pt datum =
          pt_ref - (pt_ref * normal_plane - distance_plane) * normal_plane;

      a_self->obj_ptr->setDatum(datum);
      a_self->obj_ptr->setReferenceFrame(frame);
      a_self->obj_ptr->setAlignedParaboloid(
          IRL::AlignedParaboloid({1.0e-12, -1.0e-12}));
    }
  }
}

void c_Paraboloid_setParaboloidJibben(c_Paraboloid* a_self,
                                      c_PlanarSep* a_plane, c_RectCub* a_cell,
                                      int* a_npoly, double* a_vfrac,
                                      int* a_nvert, double* a_vert_coords) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_plane != nullptr);
  assert(a_plane->obj_ptr != nullptr);
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  a_self->obj_ptr->markAsRealReconstruction();
  IRL::Plane plane = (*a_plane->obj_ptr)[0];
  IRL::Normal normal_plane = plane.normal();
  double distance_plane = plane.distance();
  if (IRL::squaredMagnitude(normal_plane) < DBL_MIN) {
    if (distance_plane < 0.0) {
      a_self->obj_ptr->markAsAlwaysBelow();
    } else {
      a_self->obj_ptr->markAsAlwaysAbove();
    }
  } else {
    IRL::RectangularCuboid cube = (*a_cell->obj_ptr);

    /* Define scale ! */
    const double scale =
        (cube.calculateSideLength(0) + cube.calculateSideLength(1) +
         cube.calculateSideLength(2)) /
        3.0;
    IRL::UnsignedIndex_t count = 0;
    // for (IRL::UnsignedIndex_t n = 0; n < *a_npoly; ++n) {
    //   for (IRL::UnsignedIndex_t m = 0; m < a_nvert[n]; ++m) {
    //     a_vert_coords[count++] /= scale;
    //   }
    // }

    /* Creating list of polygons from coordinates */
    std::vector<IRL::Polygon> polygons;
    polygons.resize(*a_npoly);
    count = 0;
    for (IRL::UnsignedIndex_t n = 0; n < polygons.size(); ++n) {
      for (IRL::UnsignedIndex_t m = 0; m < a_nvert[n]; ++m) {
        polygons[n].addVertex(IRL::Pt(a_vert_coords[3 * count + 0] / scale,
                                      a_vert_coords[3 * count + 1] / scale,
                                      a_vert_coords[3 * count + 2] / scale));
        count++;
      }
      polygons[n].calculateAndSetPlaneOfExistence();
    }

    /* Setting up reference point */
    IRL::Pt pref = polygons[0].calculateCentroid();

    /* Setting up frame of reference */
    IRL::ReferenceFrame frame;
    IRL::UnsignedIndex_t largest_dir = 0;
    if (std::fabs(normal_plane[largest_dir]) < std::fabs(normal_plane[1]))
      largest_dir = 1;
    if (std::fabs(normal_plane[largest_dir]) < std::fabs(normal_plane[2]))
      largest_dir = 2;
    if (largest_dir == 0)
      frame[0] = crossProduct(normal_plane, IRL::Normal(0.0, 1.0, 0.0));
    else if (largest_dir == 1)
      frame[0] = crossProduct(normal_plane, IRL::Normal(0.0, 0.0, 1.0));
    else
      frame[0] = crossProduct(normal_plane, IRL::Normal(1.0, 0.0, 0.0));
    frame[0].normalize();
    frame[1] = crossProduct(normal_plane, frame[0]);
    frame[2] = normal_plane;

    Eigen::MatrixXd A_mat = Eigen::MatrixXd::Zero(6, 6);
    Eigen::VectorXd b_vec = Eigen::VectorXd::Zero(6);

    for (IRL::UnsignedIndex_t n = 0; n < polygons.size(); ++n) {
      const IRL::UnsignedIndex_t shape = polygons[n].getNumberOfVertices();
      if (shape == 0) {
        continue;
      }
      // Local polygon normal and centroid
      IRL::Pt ploc = polygons[n].calculateCentroid();
      IRL::Normal nloc = polygons[n].calculateNormal();
      if (nloc * normal_plane < -0.5) {
        continue;
      }
      ploc -= pref;
      const IRL::Pt tmp_pt = ploc;
      const IRL::Normal tmp_n = nloc;
      for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
        ploc[d] = frame[d] * tmp_pt;
        nloc[d] = frame[d] * tmp_n;
      }
      // Plane coefficients
      Eigen::VectorXd reconstruction_plane_coeffs(3);
      reconstruction_plane_coeffs << -(ploc * nloc), nloc[0], nloc[1];
      reconstruction_plane_coeffs /= -IRL::safelyTiny(nloc[2]);
      // Integrals
      Eigen::VectorXd integrals = Eigen::VectorXd::Zero(6);
      double b_dot_sum = 0.0;
      for (IRL::UnsignedIndex_t v = 0; v < shape; ++v) {
        IRL::UnsignedIndex_t vn = (v + 1) % shape;
        IRL::Pt vert1 = polygons[n][v];
        IRL::Pt vert2 = polygons[n][vn];
        vert1 -= pref;
        vert2 -= pref;
        IRL::Pt tmp_pt1 = vert1;
        IRL::Pt tmp_pt2 = vert2;
        for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
          vert1[d] = frame[d] * tmp_pt1;
          vert2[d] = frame[d] * tmp_pt2;
        }

        const double xv = vert1[0];
        const double yv = vert1[1];
        const double xvn = vert2[0];
        const double yvn = vert2[1];

        Eigen::VectorXd integral_to_add(6);
        integral_to_add << (xv * yvn - xvn * yv) / 2.0,
            (xv + xvn) * (xv * yvn - xvn * yv) / 6.0,
            (yv + yvn) * (xv * yvn - xvn * yv) / 6.0,
            (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0,
            (yvn - yv) *
                (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
                 2.0 * xv * xvn * yvn + xvn * xvn * yv +
                 3.0 * xvn * xvn * yvn) /
                24.0,
            (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
        integrals += integral_to_add;
      }
      b_dot_sum += integrals.head(3).dot(reconstruction_plane_coeffs);

      // Get weighting
      const double vfrac = a_vfrac[n];
      double vfrac_weight = 1.0;
      const double limit_vfrac = 0.1;
      if (vfrac < limit_vfrac) {
        vfrac_weight = 0.5 - 0.5 * std::cos(M_PI * vfrac / limit_vfrac);
      } else if (vfrac > 1.0 - limit_vfrac) {
        vfrac_weight = 0.5 - 0.5 * std::cos(M_PI * (1.0 - vfrac) / limit_vfrac);
      }
      double ww = 1.0;
      ww *= vfrac_weight;

      if (ww > 0.0) {
        A_mat += ww * integrals * integrals.transpose();
        b_vec += ww * integrals * b_dot_sum;
      }
    }

    Eigen::VectorXd sol = A_mat.colPivHouseholderQr().solve(b_vec);
    const double a = sol(0), b = sol(1), c = sol(2), d = sol(3), e = sol(4),
                 f = sol(5);
    const double theta = 0.5 * std::atan2(e, (IRL::safelyTiny(d - f)));
    const double cos_t = std::cos(theta);
    const double sin_t = std::sin(theta);
    double A = -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
    double B = -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);

    // Translation to coordinate system R' where aligned paraboloid
    // valid Translation is R' = {x' = x + u, y' = y + v, z' = z + w}
    const double denominator = IRL::safelyTiny(4.0 * d * f - e * e);
    double u = (2.0 * b * f - c * e) / denominator;
    double v = -(b * e - 2.0 * d * c) / denominator;
    double w = -(a + (-b * b * f + b * c * e - c * c * d) / denominator);

    IRL::UnitQuaternion rotation(theta, frame[2]);
    IRL::Pt datum = pref - u * frame[0] - v * frame[1] - w * frame[2];
    auto new_frame = rotation * frame;

    if (std::isnan(u)) {
      std::cout << "u ISNAN" << std::endl;
    }
    if (std::isnan(v)) {
      std::cout << "v ISNAN" << std::endl;
    }
    if (std::isnan(w)) {
      std::cout << "v ISNAN" << std::endl;
    }
    if (std::isnan(A)) {
      std::cout << "A ISNAN" << std::endl;
    }
    if (std::isnan(A)) {
      std::cout << "B ISNAN" << std::endl;
    }

    /* Rescale everything */
    datum *= scale;
    pref *= scale;
    u *= scale;
    v *= scale;
    w *= scale;
    A /= scale;
    B /= scale;

    if (std::fabs(A) < 1.0e-6) {
      A = std::copysign(1.0e-6, A);
    }
    if (std::fabs(B) < 1.0e-6) {
      B = std::copysign(1.0e-6, B);
    }

    const double max_curvature_dx = 2.0;
    if (std::sqrt(u * u + v * v + w * w) > 20.0 * scale ||
        std::fabs(A) * scale > max_curvature_dx ||
        std::fabs(B) * scale > max_curvature_dx) {
      A = 1.0e-6;
      B = -1.0e-6;
      datum = pref - (pref * normal_plane - plane.distance()) * normal_plane;
    }

    a_self->obj_ptr->setDatum(datum);
    a_self->obj_ptr->setReferenceFrame(new_frame);
    a_self->obj_ptr->setAlignedParaboloid(IRL::AlignedParaboloid({A, B}));

    const auto cube_copy = cube;
    IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
        solver_distance(cube, a_vfrac[0], 1.0e-14,
                        IRL::Paraboloid(datum, new_frame, A, B));

    // std::cout << std::scientific << std::setprecision(15)
    //           << "Distance = " << solver_distance.getDistance() <<
    // std::endl;
    if (solver_distance.getDistance() < -999.0) {
      std::cout << "SOLVER DISTANCE FAILED" << std::endl;
      datum = pref - (pref * normal_plane - plane.distance()) * normal_plane;
      a_self->obj_ptr->setDatum(datum);
      a_self->obj_ptr->setReferenceFrame(new_frame);
      a_self->obj_ptr->setAlignedParaboloid(
          IRL::AlignedParaboloid({1.0e-12, -1.0e-12}));
    } else if (std::isnan(solver_distance.getDistance())) {
      std::cout << "SOLVER DISTANCE ISNAN" << std::endl;
    } else {
      auto new_datum =
          IRL::Pt(datum + solver_distance.getDistance() * new_frame[2]);
      a_self->obj_ptr->setDatum(new_datum);
    }

    auto volume = IRL::getVolumeMoments<IRL::Volume>(
        cube_copy, IRL::Paraboloid(*a_self->obj_ptr));

    if (std::fabs(volume / cube.calculateVolume() - a_vfrac[0]) > 1.0e-13) {
      std::cout << "Error: vfrac = " << volume / cube.calculateVolume()
                << " instead of " << a_vfrac[0] << "   ERROR = "
                << std::fabs(volume / cube.calculateVolume() - a_vfrac[0])
                << std::endl;
      std::cout << cube << std::endl;
      std::cout << (*a_self->obj_ptr) << std::endl;
      exit(1);
    }

    // if (magnitude(a_self->obj_ptr->getDatum() -
    //               IRL::Pt(-0.0001025, -0.000165, 0.0001588)) < 1.0e-3 / 25.0)
    //   std::cout << (*a_self->obj_ptr) << std::endl;
  }
}

void c_Paraboloid_setParaboloidFull(c_Paraboloid* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->markAsAlwaysAbove();
}

void c_Paraboloid_setParaboloidEmpty(c_Paraboloid* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->markAsAlwaysBelow();
}

void c_Paraboloid_copy(c_Paraboloid* a_self,
                       const c_Paraboloid* a_other_planar_separator) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_other_planar_separator != nullptr);
  assert(a_other_planar_separator->obj_ptr != nullptr);
  (*a_self->obj_ptr) = (*a_other_planar_separator->obj_ptr);
}

void c_Paraboloid_getDatum(c_Paraboloid* a_self, double* a_datum) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_datum[0] = (*a_self->obj_ptr).getDatum()[0];
  a_datum[1] = (*a_self->obj_ptr).getDatum()[1];
  a_datum[2] = (*a_self->obj_ptr).getDatum()[2];
}

void c_Paraboloid_getReferenceFrame(c_Paraboloid* a_self, double* a_frame) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_frame[0] = (*a_self->obj_ptr).getReferenceFrame()[0][0];
  a_frame[1] = (*a_self->obj_ptr).getReferenceFrame()[0][1];
  a_frame[2] = (*a_self->obj_ptr).getReferenceFrame()[0][2];
  a_frame[3] = (*a_self->obj_ptr).getReferenceFrame()[1][0];
  a_frame[4] = (*a_self->obj_ptr).getReferenceFrame()[1][1];
  a_frame[5] = (*a_self->obj_ptr).getReferenceFrame()[1][2];
  a_frame[6] = (*a_self->obj_ptr).getReferenceFrame()[2][0];
  a_frame[7] = (*a_self->obj_ptr).getReferenceFrame()[2][1];
  a_frame[8] = (*a_self->obj_ptr).getReferenceFrame()[2][2];
}

void c_Paraboloid_getAlignedParaboloid(c_Paraboloid* a_self,
                                       double* a_aligned_paraboloid) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_aligned_paraboloid[0] = (*a_self->obj_ptr).getAlignedParaboloid().a();
  a_aligned_paraboloid[1] = (*a_self->obj_ptr).getAlignedParaboloid().b();
}

void c_Paraboloid_triangulateInsideCuboid(c_Paraboloid* a_self,
                                          c_RectCub* a_cell,
                                          c_TriangulatedParaboloid* a_surface) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_cell != nullptr);
  assert(a_cell->obj_ptr != nullptr);
  assert(a_surface != nullptr);
  assert(a_surface->obj_ptr != nullptr);
  IRL::RectangularCuboid cube = (*a_cell->obj_ptr);
  IRL::Paraboloid paraboloid = (*a_self->obj_ptr);
  auto moments = IRL::getVolumeMoments<
      IRL::AddSurfaceOutput<IRL::Volume, IRL::ParametrizedSurfaceOutput>>(
      cube, paraboloid);
  double length_scale = std::cbrt(cube.calculateVolume()) / 3.0;
  moments.getSurface().triangulate_fromPtr(length_scale, 3, a_surface->obj_ptr);
  // auto surface = moments.getSurface().triangulate(length_scale, 3);
  // (*a_surface->obj_ptr).getVertexList() = std::move(surface.getVertexList());
}

void c_Paraboloid_printToScreen(const c_Paraboloid* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << (*a_self->obj_ptr);
}

}  // end extern C
