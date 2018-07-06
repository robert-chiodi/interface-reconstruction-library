/*
 * rectangular_cuboid_cut.cpp
 *
 *  Created on: Jul 5, 2018
 *      Author: Robert Chiodi
 */

#include "src/geometric_cutting.h"

namespace R2P {

RectangularCuboid getShiftedUnitCube(double a_x_shift, double a_y_shift,
                                     double a_z_shift) {
  return RectangularCuboid(
      Pt(-0.5 + a_x_shift, -0.5 + a_y_shift, -0.5 + a_z_shift),
      Pt(0.5 + a_x_shift, 0.5 + a_y_shift, 0.5 + a_z_shift));
}

void getCutTet(const Tet& a_tet, const double* a_vert_distances,
               const int a_cutting_case, Pt* a_vertices) {
  a_vertices[0] = a_tet.vertex_m[0];
  a_vertices[1] = a_tet.vertex_m[1];
  a_vertices[2] = a_tet.vertex_m[2];
  a_vertices[3] = a_tet.vertex_m[3];
  for (int v = 0; v < number_of_new_vertices_after_cut[a_cutting_case]; ++v) {
    int v1 = cut_vertices[a_cutting_case][0][v];
    int v2 = cut_vertices[a_cutting_case][1][v];
    a_vertices[4 + v] =
        pointPlaneIntersectsLine(a_vertices[v1], a_vert_distances[v1],
                                 a_vertices[v2], a_vert_distances[v2]);
  }
}

Tet tetFromRectangularCuboid(const RectangularCuboid& a_rectangular_cuboid,
                             const int a_tet_number) {
  return Tet(a_rectangular_cuboid
                 .vertex_m[rectangular_cuboid_to_5_tets[a_tet_number][0]],
             a_rectangular_cuboid
                 .vertex_m[rectangular_cuboid_to_5_tets[a_tet_number][1]],
             a_rectangular_cuboid
                 .vertex_m[rectangular_cuboid_to_5_tets[a_tet_number][2]],
             a_rectangular_cuboid
                 .vertex_m[rectangular_cuboid_to_5_tets[a_tet_number][3]]);
}

Tet tetFromCutTetVertices(const Pt* a_cut_tet_vertices,
                          const int a_cutting_case,
                          const int a_tet_number_to_get) {
  return Tet(a_cut_tet_vertices[verts_for_tets[a_cutting_case]
                                              [a_tet_number_to_get][0]],
             a_cut_tet_vertices[verts_for_tets[a_cutting_case]
                                              [a_tet_number_to_get][1]],
             a_cut_tet_vertices[verts_for_tets[a_cutting_case]
                                              [a_tet_number_to_get][2]],
             a_cut_tet_vertices[verts_for_tets[a_cutting_case]
                                              [a_tet_number_to_get][3]]);
}

double distanceBetweenPoints(const Pt& a_pt_1, const Pt& a_pt_2) {
  return std::sqrt((a_pt_1.x_m - a_pt_2.x_m) * (a_pt_1.x_m - a_pt_2.x_m) +
                   (a_pt_1.y_m - a_pt_2.y_m) * (a_pt_1.y_m - a_pt_2.y_m) +
                   (a_pt_1.z_m - a_pt_2.z_m) * (a_pt_1.z_m - a_pt_2.z_m));
}

double signedDistanceToPoint(const Pt& a_point, const Plane& a_cutting_plane) {
  return a_cutting_plane.normal_m[0] * a_point.x_m +
         a_cutting_plane.normal_m[1] * a_point.y_m +
         a_cutting_plane.normal_m[2] * a_point.z_m - a_cutting_plane.distance_m;
}

void signedDistanceToTetPoints(const Tet& a_tet, const Plane& a_cutting_plane,
                               const double flip_cut,
                               double* a_distance_to_vertices) {
  for (int v = 0; v < 4; ++v) {
    a_distance_to_vertices[v] =
        flip_cut * signedDistanceToPoint(a_tet.vertex_m[v], a_cutting_plane);
  }
}

Pt pointPlaneIntersectsLine(const Pt& a_point_1, const double a_dist_1,
                            const Pt& a_point_2, const double a_dist_2) {
  double mu = std::fmin(
      1.0,
      std::fmax(0.0, -a_dist_1 /
                         std::copysign(std::fabs(a_dist_2 - a_dist_1) + DBL_MIN,
                                       a_dist_2 - a_dist_1)));
  return (1.0 - mu) * a_point_1 + mu * a_point_2;
}

int getGeometricCaseId(const double* a_distances,
                       const int a_number_of_distances) {
  int return_case_id = 0;
  for (int d = 0; d < a_number_of_distances; ++d) {
    return_case_id +=
        static_cast<int>(0.5 + std::copysign(0.51, a_distances[d])) *
        std::pow(2, d);
  }
  return return_case_id;
}

PhaseMoments cutRectangularCuboidByPlanes(
    const RectangularCuboid a_rectangular_cuboid,
    const Reconstruction& a_reconstruction) {
  // Separate rectangular cuboid into 5 tets and accumulate moments
  GeometricMoments under_plane_moments;
  for (int t = 0; t < 5; ++t) {
    under_plane_moments += cutTetByPlanes(
        tetFromRectangularCuboid(a_rectangular_cuboid, t), a_reconstruction,
        std::max(a_reconstruction.number_of_interfaces_m, 1) - 1);
  }
  // Careful, at this point under_plane_moments.centroid = vol*centroid
  // tmp_phase_moments will have actual centroid after
  // getBothPhaseMomentsFromComplement
  PhaseMoments tmp_phase_moments;
  tmp_phase_moments = getBothPhaseMomentsFromComplement<RectangularCuboid>(
      under_plane_moments, a_rectangular_cuboid);

  return getCorrectPhaseMomentsOrdering(tmp_phase_moments,
                                        a_reconstruction.isFlipped());
}

PhaseMoments getCorrectPhaseMomentsOrdering(const PhaseMoments& a_phase_moments,
                                            const bool a_flipped) {
  return a_flipped
             ? PhaseMoments(a_phase_moments.gas_m, a_phase_moments.liquid_m)
             : a_phase_moments;
}  // namespace R2P

GeometricMoments cutTetByPlanes(const Tet& a_tet,
                                const Reconstruction& a_reconstruction,
                                const int a_cutting_plane_number) {
  double distance_to_vertices[4];
  signedDistanceToTetPoints(a_tet,
                            a_reconstruction.planes_m[a_cutting_plane_number],
                            a_reconstruction.flip_cut_m, distance_to_vertices);

  int cutting_case = getGeometricCaseId(distance_to_vertices, 4);
  Pt vertices[8];
  getCutTet(a_tet, distance_to_vertices, cutting_case, vertices);

  // Create temporary tet to perform work on and moments to return
  GeometricMoments moments_to_return;
  Tet tmp_tet;

  // If cut by last plane in reconstruction, calculate moments as sum of tets
  // that are all one phase (lying under the plane)
  if (a_cutting_plane_number == 0) {
    for (int t = number_of_tets_after_cut[cutting_case];
         t > number_of_negative_tets_after_cut[cutting_case]; --t) {
      tmp_tet = tetFromCutTetVertices(vertices, cutting_case, t);
      double vol = tmp_tet.volume();
      moments_to_return.volume_m += vol;
      moments_to_return.centroid_m += (vol * tmp_tet.centroid());
    }
    return moments_to_return;
  }

  // Otherwise, more planes to cut by
  // Cut new tets that lay below plane and recursively call this function
  for (int t = number_of_tets_after_cut[cutting_case];
       t > number_of_negative_tets_after_cut[cutting_case]; --t) {
    tmp_tet = tetFromCutTetVertices(vertices, cutting_case, t);
    if (tmp_tet.volume() < DBL_EPSILON) {
      continue;
    }
    moments_to_return +=
        cutTetByPlanes(tmp_tet, a_reconstruction, a_cutting_plane_number - 1);
  }
  return moments_to_return;
}

}  // namespace R2P
