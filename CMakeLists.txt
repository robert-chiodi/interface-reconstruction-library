# Directory to install objects
# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.11...3.15)


# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(IRL VERSION 0.1.0
	LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()

# Inform CMake where to look for CMake files.
set(IRL_INSTALL_DIR "${PROJECT_SOURCE_DIR}/install/${CMAKE_BUILD_TYPE}")
set(CMAKE_INSTALL_PREFIX "${IRL_INSTALL_DIR}/absl")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules"
		      "${PROJECT_SOURCE_DIR}/cmake/Utils")

# Check valid build type
# Validate a user-supplied value for CMAKE_BUILD_TYPE.
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UC)
if(CMAKE_BUILD_TYPE)
  if(NOT "${CMAKE_BUILD_TYPE_UC}" MATCHES "^(DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL)$")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be {Debug|Release|RelWithDebInfo|MinSizeRel}")
  endif()
endif()

# Directory to install objects
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UC)

# Find external packages needed by IRL
find_package(Eigen) # Provide -D EIGEN_PATH=/path/to/Eigen

# Add packages IRL stores in external
if(USE_ABSL MATCHES "OFF")
  set(IRL_CXX_FLAGS "${IRL_CXX_FLAGS} -D IRL_NO_ABSL")
else()
  set(IRL_USE_ABSL ON)
endif()

if(IRL_USE_ABSL)
  set(ABSEIL_DIR "${PROJECT_SOURCE_DIR}/external/abseil-cpp")
  set(CMAKE_CXX_FLAGS_ORIG "${CMAKE_CXX_FLAGS}")
  if(NOT "${CMAKE_BUILD_TYPE_UC}" MATCHES "DEBUG")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${IRL_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS  "-g -O0")
  endif()
  set(BUILD_TESTING_ORIG "${BUILD_TESTING}")
  set(BUILD_TESTING OFF)
  set(CMAKE_INSTALL_PREFIX_ORIG "${CMAKE_INSTALL_PREFIX}")
  set(CMAKE_INSTALL_PREFIX "${IRL_INSTALL_DIR}")

  add_subdirectory(${ABSEIL_DIR})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_ORIG}")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX_ORIG}")
  set(BUILD_TESTING "${BUILD_TESTING_ORIG}")
endif()

# Add IRL Libraries that we are creating
add_library(irl STATIC)
add_library(irl_c STATIC)
add_library(irl_fortran STATIC)

set_property(TARGET irl PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${IRL_INSTALL_DIR}/lib")
set_property(TARGET irl PROPERTY CXX_STANDARD 14)
set_property(TARGET irl PROPERTY COMPILE_FLAGS "${IRL_CXX_FLAGS}")
set_property(TARGET irl_c PROPERTY CXX_STANDARD 14)
set_property(TARGET irl_c PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${IRL_INSTALL_DIR}/lib")
set_property(TARGET irl_c PROPERTY COMPILE_FLAGS "${IRL_CXX_FLAGS}")

# C++ IRL Base
target_include_directories(irl PRIVATE ./)

# C Interface
target_include_directories(irl_c PRIVATE ./)
target_link_libraries(irl_c PRIVATE irl)

# Fortran Interface
set_target_properties(irl_fortran 
		      PROPERTIES Fortran_MODULE_DIRECTORY "${IRL_INSTALL_DIR}/include")
set_property(TARGET irl_fortran PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${IRL_INSTALL_DIR}/lib")
set_property(TARGET irl_fortran PROPERTY COMPILE_FLAGS "${IRL_Fortran_FLAGS}")
target_include_directories(irl_fortran PRIVATE "${IRL_INSTALL_DIR}/include")
target_link_libraries(irl_fortran PRIVATE irl)
target_link_libraries(irl_fortran PRIVATE irl_c)

# Link used portions of absl
if(IRL_USE_ABSL)
  target_include_directories(irl SYSTEM PRIVATE "${ABSEIL_DIR}")
  target_include_directories(irl_c SYSTEM PRIVATE "${ABSEIL_DIR}")
  target_link_libraries(irl PRIVATE absl::inlined_vector absl::flat_hash_map)
  target_link_libraries(irl_c PRIVATE absl::inlined_vector absl::flat_hash_map)
endif()

# Directory for test files
# Will automatically download GoogleTest and compile it to
# use for the testing framework.
if("${BUILD_TESTING}" MATCHES "ON")
  include(CTest)
  enable_testing()
  add_subdirectory(./tests)
endif()

# Need to find and link to Eigen package for C++ and C.
target_include_directories(irl SYSTEM PRIVATE "${EIGEN_INCLUDE_DIR}")
target_include_directories(irl_c SYSTEM PRIVATE "${EIGEN_INCLUDE_DIR}")

# Directory with source files.
add_subdirectory("${PROJECT_SOURCE_DIR}/src")

# Directory with examples
set(IRL_EXAMPLES_OUTPUT_DIR "${IRL_INSTALL_DIR}/examples")
add_subdirectory("${PROJECT_SOURCE_DIR}/examples")

# Directory with development tools.
#add_subdirectory("${PROJECT_SOURCE_DIR}/tools")

# Export targets to be used by other programs
install(TARGETS irl EXPORT "irl_exp" ARCHIVE DESTINATION "${IRL_INSTALL_DIR}/lib")
install(TARGETS irl_c EXPORT irl_c_exp ARCHIVE DESTINATION "${IRL_INSTALL_DIR}/lib")
install(TARGETS irl_fortran EXPORT irl_fortran_exp ARCHIVE DESTINATION "${IRL_INSTALL_DIR}/lib")
