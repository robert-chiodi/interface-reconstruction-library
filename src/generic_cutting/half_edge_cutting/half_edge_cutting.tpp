// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_TPP_
#define SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_TPP_

namespace IRL {

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
splitHalfEdgePolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                      SegmentedHalfEdgePolyhedronType* a_clipped_polytope,
                      HalfEdgePolytopeType* a_complete_polytope,
                      const Plane& a_plane,
      		      const double a_volume_tolerance){
  splitHalfEdgePolytope(a_polytope, a_clipped_polytope, a_complete_polytope, a_plane);
  
  if(std::fabs(a_polytope->calculateVolume()) < a_volume_tolerance){
    a_polytope->setNumberOfVertices(0);
    a_polytope->setNumberOfFaces(0);
  }
  if(std::fabs(a_clipped_polytope->calculateVolume()) < a_volume_tolerance){
    a_clipped_polytope->setNumberOfVertices(0);
    a_clipped_polytope->setNumberOfFaces(0);
  }    
}

template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
splitHalfEdgePolytope(SegmentedHalfEdgePolygonType* a_polytope,
                      SegmentedHalfEdgePolygonType* a_clipped_polytope,
                      HalfEdgePolytopeType* a_complete_polytope,
                      const Plane& a_plane,
		      const double a_volume_tolerance){
  splitHalfEdgePolytope(a_polytope, a_clipped_polytope, a_complete_polytope, a_plane);
  if(std::fabs(a_polytope->calculateVolume()) < a_volume_tolerance){
    a_polytope->setNumberOfVertices(0);
    a_polytope->setNumberOfFaces(0);
  }
  if(std::fabs(a_clipped_polytope->calculateVolume()) < a_volume_tolerance){
    a_clipped_polytope->setNumberOfVertices(0);
    a_clipped_polytope->setNumberOfFaces(0);
  }    
}

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
truncateHalfEdgePolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                         HalfEdgePolytopeType* a_complete_polytope,
                         const Plane& a_plane,
			 const double a_volume_tolerance){
  truncateHalfEdgePolytope(a_polytope, a_complete_polytope, a_plane);
  if(std::fabs(a_polytope->calculateVolume()) < a_volume_tolerance){
    a_polytope->setNumberOfVertices(0);
    a_polytope->setNumberOfFaces(0);
  }
}

template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
truncateHalfEdgePolytope(SegmentedHalfEdgePolygonType* a_polytope,
                         HalfEdgePolytopeType* a_complete_polytope,
                         const Plane& a_plane,
			 const double a_volume_tolerance){
  truncateHalfEdgePolytope(a_polytope, a_complete_polytope, a_plane);
  if(std::fabs(a_polytope->calculateVolume()) < a_volume_tolerance){
    a_polytope->setNumberOfVertices(0);
    a_polytope->setNumberOfFaces(0);
  }  
}



  
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
splitHalfEdgePolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                      SegmentedHalfEdgePolyhedronType* a_clipped_polytope,
                      HalfEdgePolytopeType* a_complete_polytope,
                      const Plane& a_plane) {
  using VertexType = typename HalfEdgePolytopeType::vertex_type;
  using FaceType = typename HalfEdgePolytopeType::face_type;
  using HalfEdgeType = typename HalfEdgePolytopeType::half_edge_type;

  // Create variable for clipped portion of polytope
  a_clipped_polytope->setNumberOfFaces(0);
  a_clipped_polytope->setNumberOfVertices(0);

  // Calculate distance to all vertices and marked those to be clipped (d > 0.0)
  const int early_termination_flag =
      a_polytope->calculateAndStoreDistanceToVertices(a_plane);

  if (early_termination_flag == -1) {
    return;
  } else if (early_termination_flag == 1) {
    std::swap(*a_polytope, *a_clipped_polytope);
    return;
  }

  // Clear visitation knowledge from polytope.
  const auto starting_number_of_faces = a_polytope->getNumberOfFaces();
  for(UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f){
    (*a_polytope)[f]->markAsNotVisited();
  }

  const auto starting_number_of_vertices =
      a_polytope->getNumberOfVertices();

  // Check all edges for change in isClipped. If exists, make new vertex
  // and place new half edge.
  UnsignedIndex_t number_of_vertices_unclipped = 0;     
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    if (vertex.isClipped()) {
      a_clipped_polytope->addVertex(&vertex);
      // Check all edges
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        if (current_edge->getPreviousVertex()->isNotClipped()) {
          // Found intersection, create new half edges
          subdivideEdge(current_edge, a_polytope, a_complete_polytope);
        }
        // Move to next edge
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    } else {
      a_polytope->getVertexPointer(number_of_vertices_unclipped) =
	&vertex;
      ++number_of_vertices_unclipped;
    }
  }

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();
  const auto clipped_vertices_after_intersection = a_clipped_polytope->getNumberOfVertices();

  // Now perform actual separation, which just requires relinking
  // some half edges and generating new ones for the clipped polytope.
  // Only traverse on the portion being clipped.
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    // Only want to start with ones that have previousHalfEdge under plane
    // (isNotClipped()). This will be the new vertex's half edge, because that
    // is how it was set during finding in subdivideEdge.
    auto current_half_edge = vertex.getHalfEdge()->getNextHalfEdge();
    do {
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge->getVertex()->isClipped());

    // Creating the clipped polytope follows Take - Copy - Make
    // Take all half edges except for the last one
    // Copy the last half edge (that has an intersection vertex before
    // heading back to unclipped)
    // Make a new half edge connecting the last
    // intersection vertex with the first.
    auto starting_half_edge = vertex.getHalfEdge()->getNextHalfEdge();
    auto half_edge_to_copy = current_half_edge;

    // Create new face for clipped part
    auto clipped_portion_of_face =
        a_complete_polytope->getNewFace(FaceType(starting_half_edge));
    clipped_portion_of_face->markAsVisited();
    a_clipped_polytope->addFace(clipped_portion_of_face);

    // Create second intersection vertex. Will take the first from another face
    // after this loop.
    auto second_intersection = a_complete_polytope->getNewVertex(
                               VertexType(current_half_edge->getVertex()->getLocation()));
    second_intersection->markToBeClipped();
    a_clipped_polytope->addVertex(second_intersection);

    // Copied half edge that ends on intersection vertex
    auto copied_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
        second_intersection, half_edge_to_copy->getPreviousHalfEdge(), nullptr,
        clipped_portion_of_face));
    copied_half_edge->setOppositeHalfEdge(
        half_edge_to_copy->getOppositeHalfEdge());
    copied_half_edge->getOppositeHalfEdge()->setOppositeHalfEdge(
        copied_half_edge);
    copied_half_edge->getPreviousHalfEdge()->setNextHalfEdge(copied_half_edge);

    // Set vertex half edge to the new copied_half_edge
    second_intersection->setHalfEdge(copied_half_edge);

    // New half edge making connection and closing off face
    // Will set vertex later, will be second_intersection from another face
    auto new_connection_half_edge = a_complete_polytope->getNewHalfEdge(
        HalfEdgeType(nullptr, copied_half_edge, starting_half_edge,
                     clipped_portion_of_face));
    starting_half_edge->setPreviousHalfEdge(new_connection_half_edge);

    // Set connection from copied -> new.
    copied_half_edge->setNextHalfEdge(new_connection_half_edge);

    // Reset the face pointer for all half edges in the clipped face
    auto half_edge_on_clipped_face = starting_half_edge;
    do {
      half_edge_on_clipped_face->setFace(clipped_portion_of_face);
      half_edge_on_clipped_face = half_edge_on_clipped_face->getNextHalfEdge();
    } while (half_edge_on_clipped_face != starting_half_edge);    

    // And finally close off the unclipped portion.
    doubleLinkHalfEdges(vertex.getHalfEdge(), current_half_edge);
    vertex.getHalfEdge()->getFace()->setStartingHalfEdge(vertex.getHalfEdge());

    auto new_opposite = a_complete_polytope->getNewHalfEdge(
				       HalfEdgeType(&vertex,nullptr,nullptr,nullptr));
    current_half_edge->setOppositeHalfEdge(new_opposite);
    new_opposite->setOppositeHalfEdge(current_half_edge);
    doubleLinkHalfEdges(vertex.getHalfEdge(), current_half_edge);    
    current_half_edge->getFace()->setStartingHalfEdge(current_half_edge);
    vertex.setHalfEdge(new_opposite);        
  }
  
  // Now connect newly created face edges to correct vertices.
  // (those copied for the clipped polytope)
  const auto clipped_faces = a_clipped_polytope->getNumberOfFaces();
  for(UnsignedIndex_t f = 0; f < clipped_faces; ++f){
    auto starting_half_edge = *((*a_clipped_polytope)[f]->getStartingHalfEdge());
    auto made_half_edge = starting_half_edge.getPreviousHalfEdge();    
    made_half_edge->setVertex(
        starting_half_edge.getOppositeHalfEdge()->getVertex());
    // Make half edges for face made from clipping
    auto new_opposite = a_complete_polytope->getNewHalfEdge(HalfEdgeType(made_half_edge->getPreviousVertex(),nullptr,nullptr,nullptr));
    new_opposite->setOppositeHalfEdge(made_half_edge);
    made_half_edge->setOppositeHalfEdge(new_opposite);
    new_opposite->getVertex()->setHalfEdge(new_opposite);
  }

  // Now just need to cap a_polytope and clipped_polytope.
  // Capping clipped_polytope first
  for (UnsignedIndex_t v = clipped_vertices_after_intersection;
       v < a_clipped_polytope->getNumberOfVertices();
       ++v) {
    auto current_vertex = a_clipped_polytope->getVertex(v);
    if (current_vertex->doesNotNeedToSeek()) {
      auto new_face_from_intersection = a_complete_polytope->getNewFace(FaceType(current_vertex->getHalfEdge()));
      a_clipped_polytope->addFace(new_face_from_intersection);
      new_face_from_intersection->markAsVisited();
      // current_vertex has the copied half edge for new face
      auto current_half_edge = current_vertex->getHalfEdge();
      const auto starting_half_edge = current_half_edge;
      do{
	current_half_edge->getVertex()->setToSeek();
	auto next_half_edge = current_half_edge->getOppositeHalfEdge()->getVertex()->getHalfEdge();
	current_half_edge->setFace(new_face_from_intersection);
	current_half_edge->getOppositeHalfEdge()->getFace()->markAsVisited();
	doubleLinkHalfEdges(next_half_edge, current_half_edge);
	current_half_edge = next_half_edge;
      } while(current_half_edge != starting_half_edge);
    }
  }

  // Now cap face for a_polytope (the unclipped portion).
  // We have set the new vertices to have the dangling
  // half edge, so now just need to reconnect them.
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto current_vertex = a_polytope->getVertex(v);
    if(current_vertex->doesNotNeedToSeek()){
      auto new_face_from_intersection = a_complete_polytope->getNewFace(FaceType(current_vertex->getHalfEdge()));
      a_polytope->addFace(new_face_from_intersection);

      const auto starting_vertex = current_vertex;
      do{
	current_vertex->setToSeek();
	auto current_half_edge = current_vertex->getHalfEdge();
	auto next_vertex = current_half_edge->getOppositeHalfEdge()->getVertex();
	current_half_edge->setFace(new_face_from_intersection);
	current_half_edge->getOppositeHalfEdge()->getFace()->markAsVisited();
	doubleLinkHalfEdges(next_vertex->getHalfEdge(), current_half_edge);
	current_vertex = next_vertex;
      } while(current_vertex != starting_vertex);      
    }
  }

  // Now remove vertices that are no longer used in a_polytope
  if(vertices_after_intersection > starting_number_of_vertices){  
    std::copy(&a_polytope->getVertexPointer(starting_number_of_vertices),
	      &(a_polytope->getVertexPointer(vertices_after_intersection-1))+1,
	      &a_polytope->getVertexPointer(number_of_vertices_unclipped));
  }
  a_polytope->setNumberOfVertices(number_of_vertices_unclipped+vertices_after_intersection-starting_number_of_vertices);  

  // Move fully clipped faces
  UnsignedIndex_t number_of_faces_unclipped = 0;
  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    if ((*a_polytope)[f]->hasBeenVisited()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
    } else{
      if ((*a_polytope)[f]->getStartingHalfEdge()->getVertex()->isNotClipped()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
      } else{
        a_clipped_polytope->addFace((*a_polytope)[f]);
      }
    }
  }
  if(a_polytope->getNumberOfFaces() > starting_number_of_faces){  
    std::copy(&a_polytope->getFacePointer(starting_number_of_faces),
	      &(a_polytope->getFacePointer(a_polytope->getNumberOfFaces()-1))+1,
	      &a_polytope->getFacePointer(number_of_faces_unclipped));
  }
  a_polytope->setNumberOfFaces(number_of_faces_unclipped+a_polytope->getNumberOfFaces()-starting_number_of_faces);
  assert(a_polytope->checkValidHalfEdgeStructure());
  assert(a_clipped_polytope->checkValidHalfEdgeStructure());    
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
// Splitting a half-edge polygon by a plane
////////////////////////////////////////////////////
////////////////////////////////////////////////////
// During generation of a SegmentedHalfEdgePolygon, the plane that the
// polygon exists on should have been created. Any polygon split from this
// polygon should also exist on this same plane, since that's how polygons
// work (on a 2D plane). This plane is used during the calculation of the
// moments in order to assign a consistent reference direction around which
// to determine signed normals.
template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
splitHalfEdgePolytope(SegmentedHalfEdgePolygonType* a_polytope,
                      SegmentedHalfEdgePolygonType* a_clipped_polytope,
                      HalfEdgePolytopeType* a_complete_polytope,
                      const Plane& a_plane) {
  using VertexType = typename HalfEdgePolytopeType::vertex_type;
  using FaceType = typename HalfEdgePolytopeType::face_type;
  using HalfEdgeType = typename HalfEdgePolytopeType::half_edge_type;

  // Create variable for clipped portion of polytope
  a_clipped_polytope->setNumberOfFaces(0);
  a_clipped_polytope->setNumberOfVertices(0);
  a_clipped_polytope->setPlaneOfExistence(&(a_polytope->getPlaneOfExistence()));

  // Calculate distance to all vertices and marked those to be clipped (d > 0.0)
  const int early_termination_flag =
      a_polytope->calculateAndStoreDistanceToVertices(a_plane);

  if (early_termination_flag == -1) {
    return;
  } else if (early_termination_flag == 1) {
    std::swap(*a_polytope, *a_clipped_polytope);
    return;
  }

  // Clear visitation knowledge from polytope.
  const auto starting_number_of_faces = a_polytope->getNumberOfFaces();
  for(UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f){
    (*a_polytope)[f]->markAsNotVisited();
  }

  const auto starting_number_of_vertices =
      a_polytope->getNumberOfVertices();

  // Check all edges for change in isClipped. If exists, make new vertex
  // and place new half edge.
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    if (vertex.isClipped()) {
      // Check all edges
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        if (current_edge->getPreviousVertex()->isNotClipped()) {
          // Found intersection, create new half edges
          subdivideEdge(current_edge, a_polytope, a_complete_polytope);
        }
        // Move to next edge
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    }
  }

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();

  // Now perform actual separation, which just requires relinking
  // some half edges and generating new ones for the clipped polytope.
  // Only traverse on the portion being clipped.
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    // Only want to start with ones that have previousHalfEdge under plane
    // (isNotClipped()). This will be the new vertex's half edge, because that
    // is how it was set during finding.
    auto current_half_edge = vertex.getHalfEdge()->getNextHalfEdge();
    do {
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge->getVertex()->isClipped());

    // Creating the clipped polytope follows Take - Copy - Make
    // Take all half edges except for the last one
    // Copy the last half edge (that has an intersection vertex before
    // heading back to unclipped)
    // Make a new half edge connecting the last
    // intersection vertex with the first.
    auto starting_half_edge = vertex.getHalfEdge()->getNextHalfEdge();
    auto half_edge_to_copy = current_half_edge;

    // Create new face for clipped part
    auto clipped_portion_of_face =
        a_complete_polytope->getNewFace(FaceType(starting_half_edge));
    clipped_portion_of_face->markAsVisited();
    a_clipped_polytope->addFace(clipped_portion_of_face);

    // Create second intersection vertex. Will take the first from another face
    // after this loop.
    auto second_intersection = a_complete_polytope->getNewVertex(
			       VertexType(current_half_edge->getVertex()->getLocation()));
    second_intersection->markToBeClipped();
    a_clipped_polytope->addVertex(second_intersection);

    // Copied half edge that ends on intersection vertex
    auto copied_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
        second_intersection, half_edge_to_copy->getPreviousHalfEdge(), nullptr,
        clipped_portion_of_face));
    copied_half_edge->setOppositeHalfEdge(
        half_edge_to_copy->getOppositeHalfEdge());
    copied_half_edge->getOppositeHalfEdge()->setOppositeHalfEdge(
        copied_half_edge);
    copied_half_edge->getPreviousHalfEdge()->setNextHalfEdge(copied_half_edge);

    // Set vertex half edge to the new copied_half_edge
    second_intersection->setHalfEdge(copied_half_edge);

    // New half edge making connection and closing off face
    auto new_connection_half_edge = a_complete_polytope->getNewHalfEdge(
        HalfEdgeType(nullptr, copied_half_edge, starting_half_edge,
                     clipped_portion_of_face));
    starting_half_edge->setPreviousHalfEdge(new_connection_half_edge);

    // Set connection from copied -> new.
    copied_half_edge->setNextHalfEdge(new_connection_half_edge);

    // And finally close off the unclipped portion.
    doubleLinkHalfEdges(vertex.getHalfEdge(), current_half_edge);
    vertex.getHalfEdge()->getFace()->setStartingHalfEdge(vertex.getHalfEdge());

    // Reset the face pointer for all half edges in the clipped face
    auto new_clipped_face =
        starting_half_edge->getFace() != &getOpenBoundaryFace<FaceType>()
            ? clipped_portion_of_face
            : &getOpenBoundaryFace<FaceType>();
    auto half_edge_on_clipped_face = starting_half_edge;
    do {
      half_edge_on_clipped_face->setFace(new_clipped_face);
      half_edge_on_clipped_face = half_edge_on_clipped_face->getNextHalfEdge();
    } while (half_edge_on_clipped_face != starting_half_edge);
  }

  // Now connect newly created face edges to correct vertices.
  // Note, if there are OPEN_BOUNDARY_FACEs on the polygon,
  // the polytope will not be valid until after this loop since
  // a face will point to a half edge that points to the OPEN_BOUNDARY_FACE
  for (UnsignedIndex_t f = a_clipped_polytope->getNumberOfFaces() - 1;
       f != static_cast<UnsignedIndex_t>(-1); --f) {
    auto starting_half_edge = (*a_clipped_polytope)[f]->getStartingHalfEdge();
    starting_half_edge->getPreviousHalfEdge()->setVertex(
        starting_half_edge->getOppositeHalfEdge()->getVertex());
    if ((*a_clipped_polytope)[f]->getStartingHalfEdge()->getFace() ==
        &getOpenBoundaryFace<FaceType>()) {
      a_clipped_polytope->removeFace(f);
    }
  }

  // Create new clipped outer boundary face
  for (UnsignedIndex_t v = 0; v < a_clipped_polytope->getNumberOfVertices();
       ++v) {
    auto current_vertex = a_clipped_polytope->getVertex(v);
    if (current_vertex->getHalfEdge()->getFace() ==
        &getOpenBoundaryFace<FaceType>()) {
      auto current_half_edge = current_vertex->getHalfEdge()
                                   ->getOppositeHalfEdge()
                                   ->getPreviousHalfEdge();
      current_half_edge->getFace()->markAsVisited();
      auto new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
          current_half_edge->getPreviousVertex(), current_vertex->getHalfEdge(),
          nullptr, &getOpenBoundaryFace<FaceType>()));
      current_vertex->getHalfEdge()->setNextHalfEdge(new_half_edge);
      new_half_edge->setOppositeHalfEdge(current_half_edge);
      current_half_edge->setOppositeHalfEdge(new_half_edge);
      new_half_edge->getVertex()->setToSeek();

      auto previous_new_half_edge = new_half_edge;
      while (current_half_edge->getPreviousHalfEdge()
                 ->getOppositeHalfEdge()
                 ->getFace() != &getOpenBoundaryFace<FaceType>()) {
        current_half_edge = current_half_edge->getPreviousHalfEdge()
                                ->getOppositeHalfEdge()
                                ->getPreviousHalfEdge();
        current_half_edge->getFace()->markAsVisited();
        new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
            current_half_edge->getPreviousVertex(), previous_new_half_edge,
            nullptr, &getOpenBoundaryFace<FaceType>()));
        previous_new_half_edge->setNextHalfEdge(new_half_edge);
        new_half_edge->setOppositeHalfEdge(current_half_edge);
        current_half_edge->setOppositeHalfEdge(new_half_edge);
        new_half_edge->getVertex()->setToSeek();
        previous_new_half_edge = new_half_edge;
      }
      doubleLinkHalfEdges(
          new_half_edge,
          current_half_edge->getPreviousHalfEdge()->getOppositeHalfEdge());
    }
  }

  // Now cap face for a_polytope (the unclipped portion).
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto current_vertex = a_polytope->getVertex(v);
    if (current_vertex->getHalfEdge()->getFace() ==
        &getOpenBoundaryFace<FaceType>()) {
      auto current_half_edge = current_vertex->getHalfEdge()
                                   ->getOppositeHalfEdge()
                                   ->getPreviousHalfEdge();
      current_half_edge->getFace()->markAsVisited();
      auto new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
          current_half_edge->getPreviousVertex(), current_vertex->getHalfEdge(),
          nullptr, &getOpenBoundaryFace<FaceType>()));
      current_vertex->getHalfEdge()->setNextHalfEdge(new_half_edge);
      new_half_edge->setOppositeHalfEdge(current_half_edge);
      current_half_edge->setOppositeHalfEdge(new_half_edge);
      new_half_edge->getVertex()->setToSeek();

      auto previous_new_half_edge = new_half_edge;
      while (current_half_edge->getPreviousHalfEdge()
                 ->getOppositeHalfEdge()
                 ->getFace() != &getOpenBoundaryFace<FaceType>()) {
        current_half_edge = current_half_edge->getPreviousHalfEdge()
                                ->getOppositeHalfEdge()
                                ->getPreviousHalfEdge();
        current_half_edge->getFace()->markAsVisited();
        new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
            current_half_edge->getPreviousVertex(), previous_new_half_edge,
            nullptr, &getOpenBoundaryFace<FaceType>()));
        previous_new_half_edge->setNextHalfEdge(new_half_edge);
        new_half_edge->setOppositeHalfEdge(current_half_edge);
        current_half_edge->setOppositeHalfEdge(new_half_edge);
        new_half_edge->getVertex()->setToSeek();
        previous_new_half_edge = new_half_edge;
      }
      doubleLinkHalfEdges(
          new_half_edge,
          current_half_edge->getPreviousHalfEdge()->getOppositeHalfEdge());
    }
  }

  // Now remove vertices that are no longer used in a_polytope
  UnsignedIndex_t number_of_vertices_unclipped = 0;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    if (a_polytope->getVertex(v)->isNotClipped()) {
      a_polytope->getVertexPointer(number_of_vertices_unclipped) =
          a_polytope->getVertex(v);
      ++number_of_vertices_unclipped;
    } else {
      a_clipped_polytope->addVertex(a_polytope->getVertex(v));
    }
  }
  for(UnsignedIndex_t v = starting_number_of_vertices; v < vertices_after_intersection; ++v){
      a_polytope->getVertexPointer(number_of_vertices_unclipped) =
          a_polytope->getVertex(v);
      ++number_of_vertices_unclipped;
  }
  a_polytope->setNumberOfVertices(number_of_vertices_unclipped);

  // Move fully clipped faces
  UnsignedIndex_t number_of_faces_unclipped = 0;
  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    if ((*a_polytope)[f]->hasBeenVisited()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
    } else{
      if ((*a_polytope)[f]->getStartingHalfEdge()->getVertex()->isNotClipped()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
      } else{
        a_clipped_polytope->addFace((*a_polytope)[f]);
      }
    }
  }
  for (UnsignedIndex_t f = starting_number_of_faces; f < a_polytope->getNumberOfFaces(); ++f) {
      a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
      ++number_of_faces_unclipped;
  }
  a_polytope->setNumberOfFaces(number_of_faces_unclipped);

  assert(a_clipped_polytope->checkValidHalfEdgeStructure());
  assert(a_polytope->checkValidHalfEdgeStructure());

}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Truncating instead of splitting below this ///
//////////////////////////////////////////////////
//////////////////////////////////////////////////
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
truncateHalfEdgePolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                         HalfEdgePolytopeType* a_complete_polytope,
                         const Plane& a_plane) {
  using VertexType = typename HalfEdgePolytopeType::vertex_type;
  using FaceType = typename HalfEdgePolytopeType::face_type;
  using HalfEdgeType = typename HalfEdgePolytopeType::half_edge_type;

  // Calculate distance to all vertices and marked those to be clipped (d > 0.0)
  const int early_termination_flag =
      a_polytope->calculateAndStoreDistanceToVertices(a_plane);

  if (early_termination_flag == -1) {
    return;
  } else if (early_termination_flag == 1) {
    a_polytope->setNumberOfFaces(0);
    a_polytope->setNumberOfVertices(0);
    return;
  }

  // Clear visitation knowledge from polytope.
  const auto starting_number_of_faces = a_polytope->getNumberOfFaces();
  for(UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f){
    (*a_polytope)[f]->markAsNotVisited();
  }

  const auto starting_number_of_vertices =
      a_polytope->getNumberOfVertices();

  // Check all edges for change in isClipped. If exists, make new vertex
  // and place new half edge.
  UnsignedIndex_t number_of_vertices_unclipped = 0;    
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    if (vertex.isClipped()) {      
      // Check all edges
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        if (current_edge->getPreviousVertex()->isNotClipped()) {
          // Found intersection, create new half edges
          subdivideEdge(current_edge, a_polytope, a_complete_polytope);
        }
        // Move to next edge
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    } else {
      a_polytope->getVertexPointer(number_of_vertices_unclipped) =
	&vertex;
      ++number_of_vertices_unclipped;
    }
  }

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();

  // Now perform actual separation, which just requires relinking
  // some half edges.
  // Only traverse on the portion being clipped.
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    // Only want to start with ones that have previousHalfEdge under plane
    // (isNotClipped()). This will be the new vertex's half edge, because that
    // is how it was set during finding.
    auto current_half_edge = vertex.getHalfEdge()->getNextHalfEdge();
    do {      
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge->getVertex()->isClipped());

    // Just need to close off the unclipped portion.
    // A new edge should now be formed between the first intersection we found and
    // the second intersection we found. We'll steal the first_clipped_half_edge to do this.
    // We will also generate the opposite that will be needed for the new face by stealing
    // the half edge after the first intersection vertex
    auto new_opposite = vertex.getHalfEdge()->getNextHalfEdge();
    current_half_edge->setOppositeHalfEdge(new_opposite);
    new_opposite->setVertex(&vertex);
    new_opposite->setOppositeHalfEdge(current_half_edge);
    doubleLinkHalfEdges(vertex.getHalfEdge(), current_half_edge);    
    current_half_edge->getFace()->setStartingHalfEdge(current_half_edge);
    vertex.setHalfEdge(new_opposite);    
  }

  // Now cap face for a_polytope (the unclipped portion).
  // We have set the new vertices to have the dangling
  // half edge, so now just need to reconnect them.
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto current_vertex = a_polytope->getVertex(v);
    if(current_vertex->doesNotNeedToSeek()){
      auto new_face_from_intersection = a_complete_polytope->getNewFace(FaceType(current_vertex->getHalfEdge()));
      a_polytope->addFace(new_face_from_intersection);

      const auto starting_vertex = current_vertex;
      do{
	current_vertex->setToSeek();
	auto current_half_edge = current_vertex->getHalfEdge();
	auto next_vertex = current_half_edge->getOppositeHalfEdge()->getVertex();
	current_half_edge->setFace(new_face_from_intersection);
	current_half_edge->getOppositeHalfEdge()->getFace()->markAsVisited();
	doubleLinkHalfEdges(next_vertex->getHalfEdge(), current_half_edge);
	current_vertex = next_vertex;
      } while(current_vertex != starting_vertex);      
    }
  }

  if(vertices_after_intersection > starting_number_of_vertices){
    std::copy(&a_polytope->getVertexPointer(starting_number_of_vertices),
	      &(a_polytope->getVertexPointer(vertices_after_intersection-1))+1,
	      &a_polytope->getVertexPointer(number_of_vertices_unclipped));
  }
  a_polytope->setNumberOfVertices(number_of_vertices_unclipped+vertices_after_intersection-starting_number_of_vertices);

  // Remove fully clipped faces
  UnsignedIndex_t number_of_faces_unclipped = 0;
  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    if ((*a_polytope)[f]->hasBeenVisited()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
    } else{
      if ((*a_polytope)[f]->getStartingHalfEdge()->getVertex()->isNotClipped()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
      } 
    }
  }
  if(a_polytope->getNumberOfFaces() > starting_number_of_faces){
    std::copy(&a_polytope->getFacePointer(starting_number_of_faces),
	      &(a_polytope->getFacePointer(a_polytope->getNumberOfFaces()-1))+1,
	      &a_polytope->getFacePointer(number_of_faces_unclipped));
  }
  a_polytope->setNumberOfFaces(number_of_faces_unclipped+a_polytope->getNumberOfFaces()-starting_number_of_faces);
  assert(a_polytope->checkValidHalfEdgeStructure());
}

template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
truncateHalfEdgePolytope(SegmentedHalfEdgePolygonType* a_polytope,
                         HalfEdgePolytopeType* a_complete_polytope,
                         const Plane& a_plane) {
  using VertexType = typename HalfEdgePolytopeType::vertex_type;
  using FaceType = typename HalfEdgePolytopeType::face_type;
  using HalfEdgeType = typename HalfEdgePolytopeType::half_edge_type;

  // Calculate distance to all vertices and marked those to be clipped (d > 0.0)
  const int early_termination_flag =
      a_polytope->calculateAndStoreDistanceToVertices(a_plane);

  if (early_termination_flag == -1) {
    return;
  } else if (early_termination_flag == 1) {
    a_polytope->setNumberOfFaces(0);
    a_polytope->setNumberOfVertices(0);
    return;
  }

  // Clear visitation knowledge from polytope.
  const auto starting_number_of_faces = a_polytope->getNumberOfFaces();
  for(UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f){
    (*a_polytope)[f]->markAsNotVisited();
  }

  const auto starting_number_of_vertices =
      a_polytope->getNumberOfVertices();

  // Check all edges for change in isClipped. If exists, make new vertex
  // and place new half edge.
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    if (vertex.isClipped()) {
      // Check all edges
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        if (current_edge->getPreviousVertex()->isNotClipped()) {
          // Found intersection, create new half edges
          subdivideEdge(current_edge, a_polytope, a_complete_polytope);
        }
        // Move to next edge
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    }
  }

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();

  // Now perform actual separation, which just requires relinking
  // some half edges and generating new ones for the clipped polytope.
  // Only traverse on the portion being clipped.
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    // Only want to start with ones that have previousHalfEdge under plane
    // (isNotClipped()). This will be the new vertex's half edge, because that
    // is how it was set during finding.
    auto current_half_edge = vertex.getHalfEdge()->getNextHalfEdge();
    do {
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge->getVertex()->isClipped());

    // And finally close off the unclipped portion.
    doubleLinkHalfEdges(vertex.getHalfEdge(), current_half_edge);
    vertex.getHalfEdge()->getFace()->setStartingHalfEdge(vertex.getHalfEdge());
  }

  // Now cap face for a_polytope (the unclipped portion). Creates correct
  // OPEN_BOUNDARY_FACE circuit.
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    auto current_vertex = a_polytope->getVertex(v);
    if (current_vertex->getHalfEdge()->getFace() ==
        &getOpenBoundaryFace<FaceType>()) {
      auto current_half_edge = current_vertex->getHalfEdge()
                                   ->getOppositeHalfEdge()
                                   ->getPreviousHalfEdge();
      current_half_edge->getFace()->markAsVisited();
      auto new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
          current_half_edge->getPreviousVertex(), current_vertex->getHalfEdge(),
          nullptr, &getOpenBoundaryFace<FaceType>()));
      current_vertex->getHalfEdge()->setNextHalfEdge(new_half_edge);
      new_half_edge->setOppositeHalfEdge(current_half_edge);
      current_half_edge->setOppositeHalfEdge(new_half_edge);
      new_half_edge->getVertex()->setToSeek();

      auto previous_new_half_edge = new_half_edge;
      while (current_half_edge->getPreviousHalfEdge()
                 ->getOppositeHalfEdge()
                 ->getFace() != &getOpenBoundaryFace<FaceType>()) {
        current_half_edge = current_half_edge->getPreviousHalfEdge()
                                ->getOppositeHalfEdge()
                                ->getPreviousHalfEdge();
        current_half_edge->getFace()->markAsVisited();
        new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
            current_half_edge->getPreviousVertex(), previous_new_half_edge,
            nullptr, &getOpenBoundaryFace<FaceType>()));
        previous_new_half_edge->setNextHalfEdge(new_half_edge);
        new_half_edge->setOppositeHalfEdge(current_half_edge);
        current_half_edge->setOppositeHalfEdge(new_half_edge);
        new_half_edge->getVertex()->setToSeek();
        previous_new_half_edge = new_half_edge;
      }
      doubleLinkHalfEdges(
          new_half_edge,
          current_half_edge->getPreviousHalfEdge()->getOppositeHalfEdge());
    }
  }

  // Now remove vertices that are no longer used in a_polytope
  UnsignedIndex_t number_of_vertices_unclipped = 0;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    if (a_polytope->getVertex(v)->isNotClipped()) {
      a_polytope->getVertexPointer(number_of_vertices_unclipped) =
          a_polytope->getVertex(v);
      ++number_of_vertices_unclipped;
    }
  }
  for(UnsignedIndex_t v = starting_number_of_vertices; v < vertices_after_intersection; ++v){
      a_polytope->getVertexPointer(number_of_vertices_unclipped) =
          a_polytope->getVertex(v);
      ++number_of_vertices_unclipped;
  }
  a_polytope->setNumberOfVertices(number_of_vertices_unclipped);

  // Remove fully clipped faces
  UnsignedIndex_t number_of_faces_unclipped = 0;
  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    if ((*a_polytope)[f]->hasBeenVisited()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
    } else{
      if ((*a_polytope)[f]->getStartingHalfEdge()->getVertex()->isNotClipped()) {
        a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
        ++number_of_faces_unclipped;
      } //else{
  	    //a_complete_polytope->freeFace((*a_polytope)[f]);
      //}
    }
  }
  for (UnsignedIndex_t f = starting_number_of_faces; f < a_polytope->getNumberOfFaces(); ++f) {
      a_polytope->getFacePointer(number_of_faces_unclipped) = (*a_polytope)[f];
      ++number_of_faces_unclipped;
  }
  a_polytope->setNumberOfFaces(number_of_faces_unclipped);
  assert(a_polytope->checkValidHalfEdgeStructure());
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_TPP_
