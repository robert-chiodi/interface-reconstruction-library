set(CMAKE_CXX_COMPILER mpic++)

IRLAddExample_cpp(
  TARGET_NAME
  paraboloid_advector
  HDRS
  "basic_mesh.h"
  "data.h"
  "deformation_3d.h"
  "translation_3d.h"
  "reconstruction_types.h"
  "solver.h"
  "vof_advection.h"
  "vtk.h"
  SRCS
  "main.cpp"
  "basic_mesh.cpp"
  "vof_advection.cpp"
  "solver.cpp"
  "reconstruction_types.cpp"
  "deformation_3d.cpp"
  "translation_3d.cpp"
  "vtk.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/paraboloid_advector"
  COLLECTION
  "examples_paraboloid_advector"
)
