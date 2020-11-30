// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_HELPERS_GEOMETRIC_CUTTING_HELPERS_TPP_
#define IRL_HELPERS_GEOMETRIC_CUTTING_HELPERS_TPP_

#include "irl/generic_cutting/general/encountered_id_list.h"

namespace IRL {

namespace geometric_cutting_helpers_details {

template <class ReconstructionType>
inline UnsignedIndex_t locatePt(const Pt& a_pt,
                                const ReconstructionType& a_reconstruction_link,
                                EncounteredIdList* a_id_list) {
  const auto& deciding_reconstruction =
      a_reconstruction_link.getFirstReconstruction();
  assert(a_reconstruction_link.isIdSet());
  a_id_list->addId(a_reconstruction_link.getId());

  for (UnsignedIndex_t n = 0; n < deciding_reconstruction.getNumberOfPlanes();

       ++n) {
    if (a_reconstruction_link.hasNeighbor(n) &&
        a_id_list->isIdPresent(a_reconstruction_link.getNeighbor(n).getId())) {
      continue;
    }

    if (deciding_reconstruction[n].signedDistanceToPoint(a_pt) > 0.0) {
      assert(a_reconstruction_link.hasNeighbor(n));
      return locatePt(a_pt, a_reconstruction_link.getNeighbor(n), a_id_list);
    }
  }
  // If made it this far, point is internal to this reconstruction_link node
  return a_reconstruction_link.getId();

}

}  // namespace geometric_cutting_helpers_details  

template <class DistanceListType>
LookupIndex_t getGeometricCaseId(const DistanceListType& a_distances) {
  assert(a_distances.size() > 0);
  unsigned int cutting_case_uns = 0;
  for (UnsignedIndex_t n = 0; n < a_distances.size(); ++n) {
    if (a_distances[n] > 0.0) {
      cutting_case_uns |= 1U << n;
    }
  }
  return static_cast<LookupIndex_t>(cutting_case_uns);
}

double distanceBetweenPts(const Pt& a_pt_0, const Pt& a_pt_1) {
  return magnitude(a_pt_0 - a_pt_1);
}

double squaredDistanceBetweenPts(const Pt& a_pt_0, const Pt& a_pt_1) {
  return squaredMagnitude(a_pt_0 - a_pt_1);
}

template <class ReconstructionType>
inline UnsignedIndex_t locatePt(
    const Pt& a_pt, const ReconstructionType& a_reconstruction_link) {
  const auto& deciding_reconstruction =
      a_reconstruction_link.getFirstReconstruction();
  EncounteredIdList id_list;
  assert(a_reconstruction_link.isIdSet());
  id_list.addId(a_reconstruction_link.getId());  

  for (UnsignedIndex_t n = 0; n < deciding_reconstruction.getNumberOfPlanes();
       ++n) {

    if (deciding_reconstruction[n].signedDistanceToPoint(a_pt) > 0.0) {
      assert(a_reconstruction_link.hasNeighbor(n));
      return geometric_cutting_helpers_details::locatePt(
          a_pt, a_reconstruction_link.getNeighbor(n), &id_list);
    }
  }
  // If made it this far, point is internal to this reconstruction_link node
  return a_reconstruction_link.getId();
}

void makeNormalFaceOutwardsFromVolume(const Pt& a_volume_centroid,
                                      const Pt& a_surface_centroid,
                                      Normal* a_normal) {
  if (((*a_normal) * Normal::fromPtNormalized(a_surface_centroid -
                                              a_volume_centroid)) < 0.0) {
    *a_normal = -(*a_normal);
  }
}

template <class GeometryType, class VertexDistanceList>
void signedDistanceToVertices(const GeometryType& a_geometry,
                              const Plane& a_cutting_plane,
                              VertexDistanceList* a_distance_to_vertices) {
  for (UnsignedIndex_t v = 0; v < a_geometry.getNumberOfVertices(); ++v) {
    (*a_distance_to_vertices)[v] =
        a_cutting_plane.signedDistanceToPoint(a_geometry[v]);
  }
}

template <class GeometryType>
bool isPlaneIntersectingCell(const Plane& a_plane,
                             const GeometryType& a_geometry) {
  std::array<double, GeometryType::number_of_vertices> distance_to_vertices;
  signedDistanceToVertices(a_geometry, a_plane, &distance_to_vertices);

  const auto cutting_case = getGeometricCaseId(distance_to_vertices);
  return cutting_case != 0 &&
         cutting_case != std::pow(2, a_geometry.getNumberOfVertices()) - 1;
}

template <class ReconstructionType>
bool isPtInternal(const Pt& a_pt, const ReconstructionType& a_reconstruction) {
  
  for (UnsignedIndex_t n = 0; n < a_reconstruction.getNumberOfPlanes(); ++n) {    
    // If any plane is above, point is external
    if (a_reconstruction.flip() *
            a_reconstruction[n].signedDistanceToPoint(a_pt) >
        0.0) {
      // If reconstruction is flipped, external is actually internal
      return a_reconstruction.isFlipped() ? true : false;
    }
  }
  // If made it this far, point below all planes.
  // This means it is internal, unless reconstruction is flipped.
  return a_reconstruction.isFlipped() ? false : true;
}

template <class EncompassingType, class ReconstructionType>
int checkPureObject(const EncompassingType& a_encompassing_polyhedron,
                    const ReconstructionType& a_reconstruction) {
  double min_dist = DBL_MAX;
  double max_dist = -DBL_MAX;
  for (int plane = 0; plane < a_reconstruction.getNumberOfPlanes(); ++plane) {
    for (int v = 0; v < a_encompassing_polyhedron.getNumberOfVertices(); ++v) {
      double pt_dist = a_reconstruction.flip() *
                       a_reconstruction[plane].signedDistanceToPoint(
                           a_encompassing_polyhedron[v]);
      min_dist = std::min(min_dist, pt_dist);
      max_dist = std::max(max_dist, pt_dist);
      if (min_dist * max_dist < 0.0) {
        return 2;
      }
    }
  }
  if (min_dist < 0.0) {
    return 0;  // Purely under planes
  } else {
    return 1;  // Purely over planes
  }
}

}  // namespace IRL

#endif // IRL_HELPERS_GEOMETRIC_CUTTING_HELPERS_TPP_
