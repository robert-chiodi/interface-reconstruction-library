/*
 * geometric_cutting.h
 *
 *  Created on: Jul 4, 2018
 *      Author: Robert Chiodi
 */

#ifndef SRC_GEOMETRIC_CUTTING_H_
#define SRC_GEOMETRIC_CUTTING_H_

#include <float.h>
#include <algorithm>
#include <cmath>

#include "geometry_classes.h"
#include "lookup_tables.h"

namespace R2P {

// Shifted unit cube
RectangularCuboid getShiftedUnitCube(double a_x_shift, double a_y_shift,
                                     double a_z_shift);

// Given a tet with signed distances to points,
// returns the tet with vertices at intersection points
void getCutTet(const Tet& a_tet, const double* a_vert_distances,
               const int a_cutting_case, Pt* a_vertices);

// Returns a tet from the supplied rectangular cuboid
Tet tetFromRectangularCuboid(const RectangularCuboid& a_rectangular_cuboid,
                             const int a_tet_number);

// Returns the individual tets that together form the original tet, but each new
// tet lies solely on one side of the plane
Tet tetFromCutTetVertices(const Pt* a_cut_tet_vertices,
                          const int a_cutting_case,
                          const int a_tet_number_to_get);

// Return the distance between two points
double distanceBetweenPoints(const Pt& a_pt_1, const Pt& a_pt_2);

// Get signed distance from a point from plane,
// negative is under plane, positive above
double signedDistanceToPoint(const Pt& a_point, const Plane& a_cutting_plane);

// Get signed distance from each point of tet from plane
void signedDistanceToTetPoints(const Tet& a_tet, const Plane& a_cutting_plane,
                               const double flip_cut,
                               double* a_distance_to_vertices);

// Return point where plane intersects line
Pt pointPlaneIntersectsLine(const Pt& a_point_1, const double a_dist_1,
                            const Pt& a_point_2, const double a_dist_2);

// Return case ID from collection of signed distances from points
int getGeometricCaseId(const double* a_distances,
                       const int a_number_of_distances);

// If cutting is flipped (a_flipped is true), return flipped liquid/gas phases
// in PhaseMoments
PhaseMoments getCorrectPhaseMomentsOrdering(const PhaseMoments& a_phase_moments,
                                            const bool a_flipped);

PhaseMoments cutRectangularCuboidByPlanes(
    const RectangularCuboid a_rectangular_cuboid,
    const Reconstruction& a_reconstruction);

// Cut tet by planes to return volume and centroid under the planes
GeometricMoments cutTetByPlanes(const Tet& a_tet,
                                const Reconstruction& a_reconstruction,
                                const int a_cutting_plane_number);

// Cut tet by planes to return just volume under the planes
double cutTetByPlanesForVolume(const Tet& a_tet,
                               const Reconstruction& a_reconstruction,
                               const int a_cutting_plane_number);

// With moments known for one phase, return the other
// GeometryType must be able to calculate its own volume through .volume()
// and its centroid through .centroid()
template <class GeometryType>
PhaseMoments getBothPhaseMomentsFromComplement(
    GeometricMoments a_known_moments,
    const GeometryType& a_encompassing_geometry);

//******************************************************************* //
//     Inlined function definitions placed below this.
//******************************************************************* //

//******************************************************************* //
//     Templated function definitions placed below this
//******************************************************************* //

template <class GeometryType>
PhaseMoments getBothPhaseMomentsFromComplement(
    GeometricMoments a_known_moments,
    const GeometryType& a_encompassing_geometry) {
  // Here, we are expecting a_known_moments.centroid = vol*centroid, and NOT the
  // centroid. Will return correct centroids for both phases
  PhaseMoments return_phase_moments;
  return_phase_moments.liquid_m = a_known_moments;
  double encompassing_geometry_volume = a_encompassing_geometry.volume();
  Pt encompassing_geometry_centroid = a_encompassing_geometry.centroid();

  return_phase_moments.gas_m.volume_m =
      encompassing_geometry_volume - return_phase_moments.liquid_m.volume_m;
  return_phase_moments.gas_m.centroid_m =
      encompassing_geometry_volume * encompassing_geometry_centroid -
      return_phase_moments.liquid_m.centroid_m;

  // Normalize liquid centroid if needed
  return_phase_moments.liquid_m.centroid_m =
      return_phase_moments.liquid_m.volume_m <
              DBL_EPSILON * encompassing_geometry_volume
          ? Pt(0.0, 0.0, 0.0)
          : return_phase_moments.liquid_m.centroid_m /
                return_phase_moments.liquid_m.volume_m;

  // Normalize gas centroid if needed
  return_phase_moments.liquid_m.centroid_m =
      return_phase_moments.liquid_m.volume_m <
              DBL_EPSILON * encompassing_geometry_volume
          ? Pt(0.0, 0.0, 0.0)
          : return_phase_moments.liquid_m.centroid_m /
                return_phase_moments.liquid_m.volume_m;

  return return_phase_moments;
}

}  // namespace R2P

#endif /* SRC_GEOMETRIC_CUTTING_H_ */
