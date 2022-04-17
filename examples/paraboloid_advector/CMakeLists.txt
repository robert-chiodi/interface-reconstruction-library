IRLAddExample_cpp(
  TARGET_NAME
    paraboloid_advector
  HDRS
    "basic_mesh.h"
    "data.h"
    "diagnostics.h"
    "solver.h"
  SRCS
    "basic_mesh.cpp"
    "diagnostics.cpp"
    "main.cpp"
    "solver.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/paraboloid_advector"
  COLLECTION
    "examples_paraboloid_advector"
)

configure_file("visualize.m" "${IRL_EXAMPLES_OUTPUT_DIR}/paraboloid_advector/visualize.m" COPYONLY)
