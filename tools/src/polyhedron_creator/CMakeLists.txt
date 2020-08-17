# Add polyhedron creator tool
add_executable(polyhedron_creator)
set_target_properties(polyhedron_creator
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${IRL_TOOLS_DIR}/bin"
    COMPILE_FLAGS "${IRL_CXX_FLAGS}"
    CXX_STANDARD 14)
set_property(TARGET irl_c PROPERTY CXX_STANDARD 14)
target_include_directories(polyhedron_creator PRIVATE "${PROJECT_SOURCE_DIR}")
target_link_libraries(polyhedron_creator irl)

target_sources(polyhedron_creator PRIVATE "${IRL_TOOLS_SRC_DIR}/polyhedron_creator/polyhedron_creator.h")
target_sources(polyhedron_creator PRIVATE "${IRL_TOOLS_SRC_DIR}/polyhedron_creator/polyhedron_creator.cpp")
