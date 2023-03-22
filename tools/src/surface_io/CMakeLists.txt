# Add polyhedron creator tool
add_executable(surface_io)
set_target_properties(surface_io
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${IRL_TOOLS_DIR}/bin"
    COMPILE_FLAGS "${IRL_CXX_FLAGS}"
    CXX_STANDARD 17)
set_property(TARGET irl_c PROPERTY CXX_STANDARD 17)
target_include_directories(surface_io PRIVATE "${PROJECT_SOURCE_DIR}")
target_link_libraries(surface_io irl)

target_sources(surface_io PRIVATE "${IRL_TOOLS_SRC_DIR}/surface_io/surface_io.h")
target_sources(surface_io PRIVATE "${IRL_TOOLS_SRC_DIR}/surface_io/surface_io.cpp")
