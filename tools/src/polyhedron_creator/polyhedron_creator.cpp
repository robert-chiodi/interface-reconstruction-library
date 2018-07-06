// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "tools/src/polyhedron_creator/polyhedron_creator.h"

#include <cassert>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// Trims a string to have no leading or trailing white space, and at most one
// space between characters
static inline std::string trim(std::string a_string) {
  return std::regex_replace(a_string, std::regex("^ +| +$|( ) +"), "$1");
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    simpleErrorHandler(
        "Program expects name of file containing faces represented "
        "by list of vertex indices.");
  }

  std::string user_provided_file_name(argv[1]);
  PolyhedronCreator creator;
  creator.createPolyhedronFile(user_provided_file_name);

  return 0;
}

static inline void writeLicense(std::ostream& out){
  out << "// This file is part of the Interface Reconstruction Library (IRL)" << std::endl;
  out << "// a library for interface reconstruction and computational geometry operations" << std::endl;
  out << "//" << std::endl;
  out << "// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>" << std::endl;
  out << "//" << std::endl;
  out << "// This Source Code Form is subject to the terms of the Mozilla Public" << std::endl;
  out << "// License, v. 2.0. If a copy of the MPL was not distributed with this" << std::endl;
  out << "// file, You can obtain one at https://mozilla.org/MPL/2.0/.\n"<<std::endl;
}

void PolyhedronCreator::createPolyhedronFile(
    const std::string& a_file_name_with_face_information) {
  this->processProvidedFileWithFaceInformation(
      a_file_name_with_face_information);
  this->createHeaderFile();
  this->writeHeaderFile();
  this->closeHeaderFile();
  this->createTPPFile();
  this->writeTPPFile();
  this->closeTPPFile();
}

void PolyhedronCreator::processProvidedFileWithFaceInformation(
    const std::string& a_file_name_with_face_information) {
  this->openSimplexInformationFile(a_file_name_with_face_information);
  this->getFileNameToWriteTo();
  this->getRequestedPolyhedronName();
  this->createListOfFaceIndices();
  this->closeSimplexInformationFile();
}

void PolyhedronCreator::createHeaderFile(void) {
  header_file_name_m = base_file_name_m + ".h";
  header_file_m.open(header_file_name_m);
  writeLicense(header_file_m);
}

void PolyhedronCreator::writeHeaderFile(void) {
  this->writeHeaderIncludeGuard();
  this->writeHeaderDependencies();
  this->writeStartOfNamespace();
  this->writeSpecializationClass();
  this->writeStorageClass();
  this->writePredefinedType();
  this->writeEndOfHeaderFile();
}

void PolyhedronCreator::closeHeaderFile(void) { header_file_m.close(); }

void PolyhedronCreator::openSimplexInformationFile(
    const std::string& a_file_name_with_simplex_information) {
  file_with_simplex_information_m.open(a_file_name_with_simplex_information);
}

void PolyhedronCreator::getFileNameToWriteTo(void) {
  std::string line;
  std::getline(file_with_simplex_information_m, line);
  base_file_name_m = trim(line);
}

void PolyhedronCreator::getRequestedPolyhedronName(void) {
  std::string line;
  std::getline(file_with_simplex_information_m, line);
  requested_polyhedron_name_m = trim(line);
}

void PolyhedronCreator::createListOfFaceIndices(void) {
  std::string line;
  while (std::getline(file_with_simplex_information_m, line)) {
    this->extractIndicesForAFace(line);
  }
  number_of_vertices_in_polyhedron_m =
      face_index_list_m.getNumberOfUniqueIndices();
}

void PolyhedronCreator::extractIndicesForAFace(
    const std::string& a_line_with_simplex_indices) {
  std::stringstream stream(a_line_with_simplex_indices);
  IRL::UnsignedIndex_t face_index = 0;
  int index;
  PolyhedronFaceIndexList::FaceIndexType face_index_group;
  while (stream >> index) {
    this->makeSureReadIndexIsNonNegative(index);
    face_index_group.push_back(static_cast<IRL::UnsignedIndex_t>(index));
    ++face_index;
  }
  face_index_list_m.addFace(face_index_group);
}

void PolyhedronCreator::makeSureReadIndexIsNonNegative(const int a_read_index) {
  if (a_read_index < 0) {
    simpleErrorHandler("Face indices in the file must be >=0");
  }
}

void PolyhedronCreator::closeSimplexInformationFile(void) {
  file_with_simplex_information_m.close();
}

void PolyhedronCreator::writeHeaderIncludeGuard(void) {
  std::string upper_case_header_name = base_file_name_m;
  std::transform(base_file_name_m.begin(), base_file_name_m.end(),
                 upper_case_header_name.begin(), ::toupper);
  std::string include_string =
      "SRC_GEOMETRY_POLYHEDRONS_" + upper_case_header_name + "_H_";
  header_file_m << "#ifndef " << include_string << std::endl;
  header_file_m << "#define " << include_string << '\n' << std::endl;
}

void PolyhedronCreator::writeHeaderDependencies(void) {
  header_file_m << "#include \"src/geometry/general/stored_vertex_access.h\""
                << std::endl;
  header_file_m
      << "#include \"src/geometry/half_edge_structures/half_edge_polyhedron.h\""
      << std::endl;
  header_file_m << "#include \"src/geometry/polyhedrons/general_polyhedron.h\""
                << std::endl;
  header_file_m << "#include \"src/geometry/polyhedrons/tet.h\"" << std::endl;
  header_file_m << "#include \"src/parameters/defined_types.h\"" << std::endl;
  header_file_m << std::endl;
}

void PolyhedronCreator::writeStartOfNamespace(void) {
  header_file_m << "namespace IRL{\n" << std::endl;
}

void PolyhedronCreator::writeSpecializationClass(void) {
  this->writeSpecializationClassHeader();
  this->writeFunctionDeclarationToGenerateHalfEdgeStructure();
  this->writeFunctionDeclarationToSetHalfEdgeStructure();
  this->writeFunctionDeclarationForNumberOfSimplices();
  this->writeFunctionDeclarationForSimplexIndicesFromDecomposition();
  this->writeFunctionDeclarationForSimplexDecomposition();
  this->writeEndOfSpecializationClass();
}

void PolyhedronCreator::writeSpecializationClassHeader(void) {
  header_file_m << "template <class Derived, class VertexType>" << std::endl;
  header_file_m << "class " << requested_polyhedron_name_m << "Specialization"
                << ": public GeneralPolyhedron<Derived, VertexType, ProxyTet<Derived>> {"
                << std::endl;
  header_file_m << "public: \n" << std::endl;
}

void PolyhedronCreator::writeFunctionDeclarationToGenerateHalfEdgeStructure(
    void) {
  header_file_m << "HalfEdgePolyhedron<VertexType> "
                   "generateHalfEdgeVersion(void) const; \n"
                << std::endl;
}

void PolyhedronCreator::writeFunctionDeclarationToSetHalfEdgeStructure(void) {
  header_file_m << "template<class HalfEdgePolyhedronType>" << std::endl;
  header_file_m << "void "
                   "setHalfEdgeVersion(HalfEdgePolyhedronType* "
                   "a_half_edge_version) const;"
                << "\n"
                << std::endl;
}

void PolyhedronCreator::writeFunctionDeclarationForNumberOfSimplices(void) {
  header_file_m << "static constexpr UnsignedIndex_t "
                   "getNumberOfSimplicesInDecomposition(void);"
                << "\n"
                << std::endl;
}

void PolyhedronCreator::
    writeFunctionDeclarationForSimplexIndicesFromDecomposition(void) {
  header_file_m << "static constexpr std::array<UnsignedIndex_t, 4> "
                   "getSimplexIndicesFromDecomposition(const UnsignedIndex_t "
                   "a_tet);"
                << "\n"
                << std::endl;
}

void PolyhedronCreator::writeFunctionDeclarationForSimplexDecomposition(void) {
  header_file_m
      << "ProxyTet<Derived> "
         "getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const;"
      << "\n"
      << std::endl;
}

void PolyhedronCreator::writeEndOfSpecializationClass(void) {
  header_file_m << "}; \n" << std::endl;
}

void PolyhedronCreator::writeStorageClass(void) {
  this->writeStorageClassHeader();
  this->writeStorageInheritance();
  this->writeEndOfStorageClass();
}

void PolyhedronCreator::writeStorageClassHeader(void) {
  std::string stored_name = "Stored" + requested_polyhedron_name_m;
  std::string self_name = stored_name + "<VertexType>";
  std::string storage_inher_name =
      "StoredVertexAccess<" + self_name + ",VertexType," +
      std::to_string(number_of_vertices_in_polyhedron_m) + ">";
  std::string specialization = requested_polyhedron_name_m + "Specialization<" +
                               self_name + ", VertexType>";
  header_file_m << "template <class VertexType>" << std::endl;
  header_file_m << "class " << stored_name << ": public " << storage_inher_name
                << ", public " << specialization << "{" << std::endl;
  header_file_m << "friend " << storage_inher_name << ";\n" << std::endl;
  header_file_m << "public: \n" << std::endl;
}

void PolyhedronCreator::writeStorageInheritance(void) {
  std::string stored_name = "Stored" + requested_polyhedron_name_m;
  std::string self_name = stored_name + "<VertexType>";
  std::string storage_inher_name =
      "StoredVertexAccess<" + self_name + ",VertexType," +
      std::to_string(number_of_vertices_in_polyhedron_m) + ">";
  header_file_m << "using " << storage_inher_name << "::StoredVertexAccess;\n"
                << std::endl;
  header_file_m << stored_name << "(void) = default;" << std::endl;
}

void PolyhedronCreator::writeEndOfStorageClass(void) {
  header_file_m << "}; \n" << std::endl;
}

void PolyhedronCreator::writePredefinedType(void) {
  header_file_m << std::endl;
  header_file_m << "// Predefined types " << std::endl;
  header_file_m << "using " << requested_polyhedron_name_m << "= Stored"
                << requested_polyhedron_name_m << "<Pt>; \n"
                << std::endl;
}

void PolyhedronCreator::writeEndOfHeaderFile(void) {
  this->writeEndOfNamespace();
  this->writeIncludeOfImplementation();
  this->writeEndOfIncludeGuard();
}

void PolyhedronCreator::writeEndOfNamespace(void) {
  header_file_m << "} // namespace IRL \n" << std::endl;
}

void PolyhedronCreator::writeIncludeOfImplementation(void) {
  implementation_file_name_m = base_file_name_m + ".tpp";
  header_file_m << "#include \"" << implementation_file_name_m << "\""
                << std::endl;
}

void PolyhedronCreator::writeEndOfIncludeGuard(void) {
  std::string upper_case_header_name = base_file_name_m;
  std::transform(base_file_name_m.begin(), base_file_name_m.end(),
                 upper_case_header_name.begin(), ::toupper);
  std::string include_string =
      "SRC_GEOMETRY_POLYHEDRONS_" + upper_case_header_name + "_H_";
  header_file_m << "#endif //" << include_string << std::endl;
}

void PolyhedronCreator::createTPPFile(void) {
  // Name already set in PolyhedronCreator::writeIncludeOfImplementation()
  implementation_file_m.open(implementation_file_name_m);
  writeLicense(implementation_file_m);
}

void PolyhedronCreator::closeTPPFile(void) { implementation_file_m.close(); }

void PolyhedronCreator::writeTPPFile(void) {
  this->writeTPPIncludeGuard();
  this->writeTPPIncludes();
  this->openTPPNamespace();
  this->writeDecomposition();
  this->writeGenerateHalfEdgeDefinition();
  this->writeSetHalfEdgeDefinition();
  this->writeNumberOfSimplicesDefinition();
  this->writeSimplexIndicesFromDecompositionDefinition();
  this->writeSimplexDecompositionDefinition();
  this->closeTPPNamespace();
  this->writeEndOfTPPIncludeGuard();
}

void PolyhedronCreator::writeTPPIncludeGuard(void) {
  std::string upper_case_header_name = base_file_name_m;
  std::transform(base_file_name_m.begin(), base_file_name_m.end(),
                 upper_case_header_name.begin(), ::toupper);
  std::string include_string =
      "SRC_GEOMETRY_POLYHEDRONS_" + upper_case_header_name + "_TPP_";
  implementation_file_m << "#ifndef " << include_string << std::endl;
  implementation_file_m << "#define " << include_string << '\n' << std::endl;
}

void PolyhedronCreator::writeTPPIncludes(void) {
  implementation_file_m << "#include <cassert>\n" << std::endl;
  implementation_file_m
      << "#include "
         "\"src/geometry/general/moment_calculation_through_simplices.h\""
      << std::endl;
  implementation_file_m
      << "#include \"src/geometry/half_edge_structures/half_edge.h\""
      << std::endl;
  implementation_file_m
      << "#include \"src/geometry/half_edge_structures/half_edge.h\""
      << std::endl;
}

void PolyhedronCreator::openTPPNamespace(void) {
  implementation_file_m << "\n namespace IRL {\n" << std::endl;
}

static void writeOutTetDecompositionArray(
    std::ofstream& a_file,
    const std::vector<std::array<IRL::UnsignedIndex_t, 3>>& a_tet_faces) {
  a_file << "static constexpr std::array<std::array<UnsignedIndex_t, 3>, "
         << a_tet_faces.size() << "> face_triangle_decomposition{{"
         << std::endl;
  for (std::size_t tri = 0; tri < a_tet_faces.size() - 1; ++tri) {
    a_file << "{";
    for (std::size_t v = 0; v < a_tet_faces[tri].size() - 1; ++v) {
      a_file << a_tet_faces[tri][v] << ", ";
    }
    a_file << a_tet_faces[tri].back() << "}, " << std::endl;
  }

  a_file << "{";
  for (std::size_t v = 0; v < a_tet_faces.back().size() - 1; ++v) {
    a_file << a_tet_faces.back()[v] << ", ";
  }
  a_file << a_tet_faces.back().back() << "}}}; " << std::endl;
}

void PolyhedronCreator::writeDecomposition(void) {
  implementation_file_m << "namespace " << base_file_name_m
                        << "_triangulation {" << std::endl;
  auto datum_index = this->findMostCommonVertex();
  auto decomp = this->setUpTetDecomposition(datum_index);
  implementation_file_m << "static constexpr UnsignedIndex_t datum_index = "
                        << datum_index << ";" << std::endl;
  writeOutTetDecompositionArray(implementation_file_m, decomp);
  implementation_file_m << "} // namespace " << base_file_name_m
                        << "_triangulation \n"
                        << std::endl;
}

IRL::UnsignedIndex_t PolyhedronCreator::findMostCommonVertex(void) {
  std::vector<IRL::UnsignedIndex_t> index_occurence(
      number_of_vertices_in_polyhedron_m);
  std::fill(index_occurence.begin(), index_occurence.end(), 0);
  for (IRL::UnsignedIndex_t n = 0; n < face_index_list_m.size(); ++n) {
    for (IRL::UnsignedIndex_t v = 0; v < face_index_list_m[n].size(); ++v) {
      ++(index_occurence[face_index_list_m[n][v]]);
    }
  }

  return static_cast<IRL::UnsignedIndex_t>(std::distance(
      index_occurence.begin(),
      std::max_element(index_occurence.begin(), index_occurence.end())));
}

std::vector<std::array<IRL::UnsignedIndex_t, 3>>
PolyhedronCreator::setUpTetDecomposition(
    const IRL::UnsignedIndex_t a_datum_index) {
  std::vector<std::array<IRL::UnsignedIndex_t, 3>> face_triangles;
  for (IRL::UnsignedIndex_t n = 0; n < face_index_list_m.size(); ++n) {
    // Don't need tet for face that has datum index
    bool has_datum_index = false;
    for (IRL::UnsignedIndex_t v = 0; v < face_index_list_m[n].size(); ++v) {
      if (face_index_list_m[n][v] == a_datum_index) {
        has_datum_index = true;
        break;
      }
    }
    if (has_datum_index) {
      continue;
    }
    for (IRL::UnsignedIndex_t v = 2; v < face_index_list_m[n].size(); ++v) {
      std::array<IRL::UnsignedIndex_t, 3> triangle{{face_index_list_m[n][0],
                                                    face_index_list_m[n][v - 1],
                                                    face_index_list_m[n][v]}};
      face_triangles.push_back(triangle);
    }
    assert(face_index_list_m[n].size() >= 3);
  }
  return face_triangles;
}

void PolyhedronCreator::functionDefinitionHeaderWriter(
    std::ofstream& a_file, const std::string& return_type,
    const std::string& function_name, const std::string& function_arguments,
    const std::string& after_function_markings) {
  a_file << "template<class Derived, class VertexType>" << std::endl;
  a_file << return_type << " ";
  a_file << requested_polyhedron_name_m +
                "Specialization<Derived, VertexType>::";
  a_file << function_name << "(";
  a_file << function_arguments << ")";
  a_file << " " << after_function_markings << "{" << std::endl;
}

void PolyhedronCreator::writeGenerateHalfEdgeDefinition(void) {
  this->functionDefinitionHeaderWriter(
      implementation_file_m, "HalfEdgePolyhedron<VertexType>",
      "generateHalfEdgeVersion", "void", "const");
  implementation_file_m << "HalfEdgePolyhedron<VertexType> half_edge_version;"
                        << std::endl;
  implementation_file_m << "this->setHalfEdgeVersion(&half_edge_version);"
                        << std::endl;
  implementation_file_m << "return half_edge_version;" << std::endl;
  implementation_file_m << "}\n" << std::endl;
}

void PolyhedronCreator::writeSetHalfEdgeDefinition(void) {
  this->functionDefinitionHeaderWriter(
      implementation_file_m, "template<class HalfEdgePolyhedronType> void",
      "setHalfEdgeVersion", "HalfEdgePolyhedronType* a_half_edge_version",
      "const");
  this->writeFunctionToSetHalfEdgeStructure();
  implementation_file_m << "}\n" << std::endl;
}

void PolyhedronCreator::writeNumberOfSimplicesDefinition(void) {
  this->functionDefinitionHeaderWriter(
      implementation_file_m, "constexpr UnsignedIndex_t",
      "getNumberOfSimplicesInDecomposition", "void", " ");
  implementation_file_m << "return static_cast<UnsignedIndex_t>(" << base_file_name_m
                        << "_triangulation::face_triangle_decomposition.size());"
                        << std::endl;
  implementation_file_m << "}\n" << std::endl;
}

void PolyhedronCreator::writeSimplexIndicesFromDecompositionDefinition(void) {
  this->functionDefinitionHeaderWriter(
      implementation_file_m, "constexpr std::array<UnsignedIndex_t, 4>",
      "getSimplexIndicesFromDecomposition", "const UnsignedIndex_t a_tet", "");
  std::string class_name = requested_polyhedron_name_m + "Specialization::";
  implementation_file_m << "assert(a_tet < " << class_name
                        << "getNumberOfSimplicesInDecomposition());"
                        << std::endl;
  implementation_file_m << "return {";
  implementation_file_m
      << base_file_name_m
      << "_triangulation::face_triangle_decomposition[a_tet][0]," << std::endl;
  implementation_file_m
      << base_file_name_m
      << "_triangulation::face_triangle_decomposition[a_tet][1]," << std::endl;
  implementation_file_m
      << base_file_name_m
      << "_triangulation::face_triangle_decomposition[a_tet][2]," << std::endl;
  implementation_file_m << base_file_name_m << "_triangulation::datum_index};"
                        << std::endl;
  implementation_file_m << "}\n" << std::endl;
}

void PolyhedronCreator::writeSimplexDecompositionDefinition(void) {
  this->functionDefinitionHeaderWriter(
      implementation_file_m, "ProxyTet<Derived>", "getSimplexFromDecomposition",
      "const UnsignedIndex_t a_tet", "const");
  implementation_file_m
      << "assert(a_tet < this->getNumberOfSimplicesInDecomposition());"
      << std::endl;

  implementation_file_m
      << "return { static_cast<const "
         "Derived&>(*this),this->getSimplexIndicesFromDecomposition(a_tet)};"
      << std::endl;
  implementation_file_m << "}\n" << std::endl;
}

void PolyhedronCreator::closeTPPNamespace(void) {
  implementation_file_m << "\n } // namespace IRL" << std::endl;
}

void PolyhedronCreator::writeFunctionToSetHalfEdgeStructure(void) {
  EdgeMappingType edge_to_index_map = this->generateMapOfEdgesToLinearIndex();
  implementation_file_m << "using HalfEdgeType = typename "
                           "HalfEdgePolyhedronType::half_edge_type;\n"
                        << std::endl;
  this->writeInitializationOfHalfEdgePolyhedron();
  this->setVertices();
  this->setFaceHalfEdges(edge_to_index_map);
}

PolyhedronCreator::EdgeMappingType
PolyhedronCreator::generateMapOfEdgesToLinearIndex(void) {
  EdgeMappingType edge_to_index_map;
  edge_to_index_map.resize(number_of_vertices_in_polyhedron_m);
  for (auto& vector : edge_to_index_map) {
    vector.resize(number_of_vertices_in_polyhedron_m);
  }
  IRL::UnsignedIndex_t current_half_edge_index = 0;
  for (IRL::UnsignedIndex_t face_index = 0;
       face_index < face_index_list_m.size(); ++face_index) {
    const auto& face = face_index_list_m[face_index];
    for (IRL::UnsignedIndex_t vertex_index = 0; vertex_index < face.size();
         ++vertex_index) {
      IRL::UnsignedIndex_t starting_vertex = face[vertex_index];
      IRL::UnsignedIndex_t ending_vertex =
          face[(vertex_index + 1) % face.size()];
      edge_to_index_map[starting_vertex][ending_vertex] =
          current_half_edge_index;
      ++current_half_edge_index;
    }
  }

  return edge_to_index_map;
}

void PolyhedronCreator::writeInitializationOfHalfEdgePolyhedron(void) {
  number_of_half_edges_m = 0;
  for (IRL::UnsignedIndex_t face_index = 0;
       face_index < face_index_list_m.size(); ++face_index) {
    const auto& face = face_index_list_m[face_index];
    number_of_half_edges_m += static_cast<IRL::UnsignedIndex_t>(face.size());
  }
  implementation_file_m << "a_half_edge_version->resize(";
  implementation_file_m << number_of_half_edges_m << ", "
                        << number_of_vertices_in_polyhedron_m << ", "
                        << face_index_list_m.size() << ");\n"
                        << std::endl;
}

void PolyhedronCreator::setVertices(void) {
  implementation_file_m << "for (UnsignedIndex_t v = 0; v < "
                        << number_of_vertices_in_polyhedron_m << "; ++v){"
                        << std::endl;
  implementation_file_m
      << "a_half_edge_version->getVertexLocation(v).setLocation((*"
         "this)[v]);"
      << std::endl;
  implementation_file_m
      << "a_half_edge_version->getVertex(v).setVertexLocation(&a_"
         "half_edge_version->getVertexLocation(v));"
      << std::endl;
  implementation_file_m << "}\n" << std::endl;
}

static IRL::UnsignedIndex_t negativeIndexWrapper(const int a_index,
                                                 const std::size_t max_val) {
  return a_index < 0 ? static_cast<IRL::UnsignedIndex_t>(
                           static_cast<int>(max_val) + a_index)
                     : static_cast<IRL::UnsignedIndex_t>(a_index);
}

static void staticArrayWriter(
    std::ofstream& a_file, const IRL::UnsignedIndex_t a_size,
    const std::string& a_array_name,
    const std::vector<IRL::UnsignedIndex_t> a_mapping) {
  a_file << "static constexpr std::array<UnsignedIndex_t, " << a_size << "> "
         << a_array_name << "{{";
  for (std::size_t n = 0; n < a_mapping.size() - 1; ++n) {
    a_file << a_mapping[n] << ", ";
  }
  a_file << a_mapping.back() << "}};" << std::endl;
}

void PolyhedronCreator::setFaceHalfEdges(
    const EdgeMappingType& a_edge_to_index_map) {
  std::string base = "a_half_edge_version";
  std::vector<IRL::UnsignedIndex_t> ending_vertex_mapping(
      number_of_half_edges_m);
  std::vector<IRL::UnsignedIndex_t> previous_half_edge_mapping(
      number_of_half_edges_m);
  std::vector<IRL::UnsignedIndex_t> next_half_edge_mapping(
      number_of_half_edges_m);
  std::vector<IRL::UnsignedIndex_t> face_mapping(number_of_half_edges_m);
  std::vector<IRL::UnsignedIndex_t> opposite_half_edge_mapping(
      number_of_half_edges_m);

  for (IRL::UnsignedIndex_t face_index = 0;
       face_index < face_index_list_m.size(); ++face_index) {
    const auto& face = face_index_list_m[face_index];
    for (IRL::UnsignedIndex_t vertex_index = 0; vertex_index < face.size();
         ++vertex_index) {
      // Get correct linear indices
      IRL::UnsignedIndex_t starting_vertex = face[vertex_index];
      IRL::UnsignedIndex_t ending_vertex =
          face[(vertex_index + 1) % face.size()];
      IRL::UnsignedIndex_t current_half_edge =
          a_edge_to_index_map[starting_vertex][ending_vertex];
      IRL::UnsignedIndex_t next_half_edge =
          a_edge_to_index_map[ending_vertex]
                             [face[(vertex_index + 2) % face.size()]];
      IRL::UnsignedIndex_t previous_half_edge =
          a_edge_to_index_map[face[negativeIndexWrapper(
              static_cast<int>(vertex_index) - 1, face.size())]]
                             [starting_vertex];
      IRL::UnsignedIndex_t opposite_half_edge =
          a_edge_to_index_map[ending_vertex][starting_vertex];

      // Add entry for mappings
      ending_vertex_mapping[current_half_edge] = ending_vertex;
      previous_half_edge_mapping[current_half_edge] = previous_half_edge;
      next_half_edge_mapping[current_half_edge] = next_half_edge;
      face_mapping[current_half_edge] = face_index;
      opposite_half_edge_mapping[current_half_edge] = opposite_half_edge;
    }
  }
  staticArrayWriter(implementation_file_m, number_of_half_edges_m,
                    "ending_vertex_mapping", ending_vertex_mapping);
  staticArrayWriter(implementation_file_m, number_of_half_edges_m,
                    "previous_half_edge_mapping", previous_half_edge_mapping);
  staticArrayWriter(implementation_file_m, number_of_half_edges_m,
                    "next_half_edge_mapping", next_half_edge_mapping);
  staticArrayWriter(implementation_file_m, number_of_half_edges_m,
                    "face_mapping", face_mapping);
  staticArrayWriter(implementation_file_m, number_of_half_edges_m,
                    "opposite_half_edge_mapping", opposite_half_edge_mapping);
  implementation_file_m
      << "for (UnsignedIndex_t n = 0; n < "
         "static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n){"
      << std::endl;
  implementation_file_m << "HalfEdgeType& current_half_edge = " << base
                        << "->getHalfEdge(n);" << std::endl;
  implementation_file_m << "current_half_edge = HalfEdgeType(";
  implementation_file_m << "&" << base
                        << "->getVertex(ending_vertex_mapping[n]), ";
  implementation_file_m << "&" << base
                        << "->getHalfEdge(previous_half_edge_mapping[n]), ";
  implementation_file_m << "&" << base
                        << "->getHalfEdge(next_half_edge_mapping[n]), ";
  implementation_file_m << "&" << base << "->getFace(face_mapping[n]));"
                        << std::endl;
  implementation_file_m << "current_half_edge.setOppositeHalfEdge(&" << base
                        << "->getHalfEdge(opposite_half_edge_mapping[n]));"
                        << std::endl;
  implementation_file_m << "current_half_edge.getFace()->setStartingHalfEdge(&"
                        << "current_half_edge);" << std::endl;
  implementation_file_m << "current_half_edge.getVertex()->setHalfEdge(&"
                        << "current_half_edge);" << std::endl;
  implementation_file_m << "}\n" << std::endl;
}

void PolyhedronCreator::writeEndOfTPPIncludeGuard(void) {
  std::string upper_case_header_name = base_file_name_m;
  std::transform(base_file_name_m.begin(), base_file_name_m.end(),
                 upper_case_header_name.begin(), ::toupper);
  std::string include_string =
      "SRC_GEOMETRY_POLYHEDRONS_" + upper_case_header_name + "_TPP_";
  implementation_file_m << "#endif //" << include_string << std::endl;
}

