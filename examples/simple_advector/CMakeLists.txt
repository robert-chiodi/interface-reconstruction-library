IRLAddExample_cpp(
  TARGET_NAME
    simple_advector
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
    "${IRL_EXAMPLES_OUTPUT_DIR}/simple_advector"
  COLLECTION
    "examples_simple_advector"
)

configure_file("visualize.m" "${IRL_EXAMPLES_OUTPUT_DIR}/simple_advector/visualize.m" COPYONLY)
