set(IRL_TOOLS_SRC_DIR "${IRL_TOOLS_DIR}/src")
# Add tool directories, which will have CMakeLists.txt
# to handle their compilation/installation.
add_subdirectory("${IRL_TOOLS_SRC_DIR}/polyhedron_creator")
add_subdirectory("${IRL_TOOLS_SRC_DIR}/cmake_list_generator")

