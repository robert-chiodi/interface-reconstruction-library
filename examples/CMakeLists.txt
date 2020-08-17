
add_custom_target(examples_cpp)
add_custom_target(examples_fortran)
add_custom_target(examples_simple_advector)
add_custom_target(examples_advector)

include(${PROJECT_SOURCE_DIR}/cmake/Utils/IRLAddExample.cmake)

add_subdirectory("./fortran")
add_subdirectory("./cpp")
add_subdirectory("./simple_advector")
add_subdirectory("./advector")
