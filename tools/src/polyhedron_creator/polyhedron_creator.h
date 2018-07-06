// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef TOOLS_SRC_POLYHEDRON_CREATOR_POLYHEDRON_CREATOR_H_
#define TOOLS_SRC_POLYHEDRON_CREATOR_POLYHEDRON_CREATOR_H_

// This file stores the classes need to read in a file containing
// information on a polyhedron decomposed into simplices and then
// write out a class that can work inside IRL.
//
// The file containing the simplex information must be formatted as followed
// name_of_file (newline)
// NameOfPolyhedronClass (newline)
// List of vertex indices making up the faces.
//
// To adhere with conventions used in IRL, the
// file name should be lower-case and separated by
// underscores, such as rectangular_cuboid.
// The polyhedron class name should have no spaces
// or underscores and start each new word with a capital letter,
// such as RectangularCuboid.
//
// An example of such an input file would be this, used to generate
// the RectangularCuboid class.
//
// rectangular_cuboid
// RectangularCuboid
// 4 7 6 5
// 0 1 2 3
// 4 0 3 7
// 2 1 5 6
// 1 0 4 5
// 3 2 6 7

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "src/parameters/defined_types.h"

int main(int argc, char* argv[]);

void simpleErrorHandler(const std::string& error_message);

class PolyhedronFaceIndexList {
 public:
  using FaceIndexType = std::vector<IRL::UnsignedIndex_t>;

  PolyhedronFaceIndexList(void) = default;

  void addFace(const FaceIndexType& a_index_group_for_face) {
    vector_of_face_indices_m.push_back(a_index_group_for_face);
  }

  const FaceIndexType& operator[](const IRL::UnsignedIndex_t a_index) const {
    return vector_of_face_indices_m[a_index];
  }

  IRL::UnsignedIndex_t size(void) const {
    return static_cast<IRL::UnsignedIndex_t>(vector_of_face_indices_m.size());
  }

  IRL::UnsignedIndex_t getNumberOfUniqueIndices(void) const {
    std::vector<IRL::UnsignedIndex_t> unique_indices;
    for (const auto& index_group : vector_of_face_indices_m) {
      for (const auto& index : index_group) {
        auto location =
            std::find(unique_indices.begin(), unique_indices.end(), index);
        if (location == unique_indices.end()) {
          // Index hasn't been encountered yet.
          unique_indices.push_back(index);
        }
      }
    }
    this->checkNoHolesInVertexNumbers(&unique_indices);
    return static_cast<IRL::UnsignedIndex_t>(unique_indices.size());
  }

  void checkNoHolesInVertexNumbers(
      std::vector<IRL::UnsignedIndex_t>* a_vertex_list) const {
    std::sort(a_vertex_list->begin(), a_vertex_list->end());
    for (IRL::UnsignedIndex_t n = 0; n < a_vertex_list->size(); ++n) {
      if (n != (*a_vertex_list)[n]) {
        simpleErrorHandler(
            "Vertices given did not contain 0 or not all vertices were used. "
            "\n\n");
      }
    }
  }

  ~PolyhedronFaceIndexList(void) = default;

 private:
  std::vector<FaceIndexType> vector_of_face_indices_m;
};

class PolyhedronCreator {
  std::string base_file_name_m;
  std::string requested_polyhedron_name_m;
  std::string header_file_name_m;
  std::string implementation_file_name_m;
  std::ifstream file_with_simplex_information_m;
  PolyhedronFaceIndexList face_index_list_m;
  IRL::UnsignedIndex_t number_of_vertices_in_polyhedron_m;
  IRL::UnsignedIndex_t number_of_half_edges_m;
  std::ofstream header_file_m;
  std::ofstream implementation_file_m;

 public:
  PolyhedronCreator(void) = default;

  void createPolyhedronFile(
      const std::string& a_file_name_with_simplex_information);

  ~PolyhedronCreator(void) = default;

 private:
  void processProvidedFileWithFaceInformation(
      const std::string& a_file_name_with_face_information);
  void openSimplexInformationFile(
      const std::string& a_file_name_with_face_information);
  void getFileNameToWriteTo();
  void getRequestedPolyhedronName(void);

  void createListOfFaceIndices(void);
  void extractIndicesForAFace(const std::string& a_line_with_simplex_indices);
  void makeSureReadIndexIsNonNegative(const int a_read_index);
  void makeSureOnlyFourIndicesReceived(
      const IRL::UnsignedIndex_t a_length_read_in);

  void closeSimplexInformationFile(void);

  void createHeaderFile(void);

  void writeHeaderFile(void);

  void writeHeaderIncludeGuard(void);
  void writeHeaderDependencies(void);
  void writeStartOfNamespace(void);

  void writeSpecializationClass(void);
  void writeSpecializationClassHeader(void);
  void writeFunctionDeclarationToGenerateHalfEdgeStructure(void);
  void writeFunctionDeclarationToSetHalfEdgeStructure(void);
  void writeFunctionDeclarationForNumberOfSimplices(void);
  void writeFunctionDeclarationForSimplexIndicesFromDecomposition(void);
  void writeFunctionDeclarationForSimplexDecomposition(void);
  void writeEndOfSpecializationClass(void);

  void writeStorageClass(void);
  void writeStorageClassHeader(void);
  void writeStorageInheritance(void);
  void writeEndOfStorageClass(void);

  void writePredefinedType(void);

  void writeEndOfHeaderFile(void);
  void writeEndOfNamespace(void);
  void writeIncludeOfImplementation(void);
  void writeEndOfIncludeGuard(void);

  void closeHeaderFile(void);

  void createTPPFile(void);

  void writeTPPFile(void);
  void writeTPPIncludeGuard(void);
  void writeTPPIncludes(void);
  void openTPPNamespace(void);
  void writeDecomposition(void);
  void writeGenerateHalfEdgeDefinition(void);
  void writeSetHalfEdgeDefinition(void);
  void writeNumberOfSimplicesDefinition(void);
  void writeSimplexIndicesFromDecompositionDefinition(void);
  void writeSimplexDecompositionDefinition(void);
  void closeTPPNamespace(void);
  void writeEndOfTPPIncludeGuard(void);

  void writeFunctionToSetHalfEdgeStructure(void);
  using EdgeMappingType = std::vector<std::vector<IRL::UnsignedIndex_t>>;
  EdgeMappingType generateMapOfEdgesToLinearIndex(void);
  void writeHalfEdgeSetFunctionDefintion(void);
  void writeInitializationOfHalfEdgePolyhedron(void);
  void setVertices(void);
  void setFaceHalfEdges(const EdgeMappingType& a_edge_to_index_map);
  void writeEndToHalfEdgeSetFunctionDefinition(void);

  IRL::UnsignedIndex_t findMostCommonVertex(void);
  std::vector<std::array<IRL::UnsignedIndex_t, 3>> setUpTetDecomposition(
      const IRL::UnsignedIndex_t a_datum_index);
  void functionDefinitionHeaderWriter(
      std::ofstream& a_file, const std::string& return_type,
      const std::string& function_name, const std::string& function_arguments,
      const std::string& after_function_markings);

  void closeTPPFile(void);
};

void simpleErrorHandler(const std::string& error_message) {
  std::cout << error_message << std::endl;
  exit(1);
}

#endif  // TOOLS_SRC_POLYHEDRON_CREATOR_POLYHEDRON_CREATOR_H_
