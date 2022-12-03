cmake_minimum_required(VERSION 3.9...3.18)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()


# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(IRL VERSION 0.1.0
	DESCRIPTION "Interface Reconstruction Library"
	LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Inform CMake where to look for CMake files.
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules"
		      "${PROJECT_SOURCE_DIR}/cmake/Utils")

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()

# Set a default build type if none was specified
set(default_build_type "Release")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Check valid build type
# Validate a user-supplied value for CMAKE_BUILD_TYPE.
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UC)
if(CMAKE_BUILD_TYPE)
  if(NOT "${CMAKE_BUILD_TYPE_UC}" MATCHES "^(DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL)$")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be {Debug|Release|RelWithDebInfo|MinSizeRel}")
  endif()
endif()

#If not debug, add -DNDEBUG and -DNDEBUG_PERF to compiler options
if(NOT "${CMAKE_BUILD_TYPE_UC}" MATCHES "DEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -DNDEBUG_PERF") 
endif()

# Find external packages needed by IRL
find_package(Eigen3) # Provide -D EIGEN_PATH=/path/to/Eigen

# Add packages IRL stores in external
set(IRL_USE_ABSL ON)
if(DEFINED USE_ABSL)
  if(NOT USE_ABSL)
    set(IRL_CXX_FLAGS "${IRL_CXX_FLAGS} -D IRL_NO_ABSL")
  endif()
endif()

if(IRL_USE_ABSL)
  # Add Abseil as a sub-directory but temporarily turn testing off.
  set(ABSL_PROPAGATE_CXX_STD ON)
  set(ABSEIL_DIR "${PROJECT_SOURCE_DIR}/external/abseil-cpp")
  set(BUILD_TESTING_ORIG ${BUILD_TESTING})
  set(BUILD_TESTING OFF)
  set(CMAKE_INSTALL_PREFIX_ORIG "${CMAKE_INSTALL_PREFIX}")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/absl")
  add_subdirectory(${ABSEIL_DIR})
  set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX_ORIG}")
  set(BUILD_TESTING ${BUILD_TESTING_ORIG})
endif()

# Path for Fortran Module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod )

# Add IRL Libraries that we are creating
if(BUILD_SHARED_LIBS)
  add_library(irl SHARED)
  add_library(irl_c SHARED)
  add_library(irl_fortran SHARED)
  if(APPLE AND IRL_USE_ABSL)
    target_link_libraries(irl PUBLIC "-framework CoreFoundation")
  endif()
else()
  add_library(irl STATIC)
  add_library(irl_c STATIC)
  add_library(irl_fortran STATIC)
endif()

# C++ IRL Base
target_include_directories(irl PUBLIC "${PROJECT_SOURCE_DIR}")
target_link_libraries(irl INTERFACE PUBLIC Eigen3::Eigen)
if(IRL_USE_ABSL)
  target_include_directories(irl SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/abseil-cpp")
  target_link_libraries(irl PUBLIC absl::inlined_vector absl::flat_hash_map)
endif()

# C Interface
target_link_libraries(irl_c PUBLIC irl)

# Fortran Interface
target_link_libraries(irl_fortran PUBLIC irl_c)

if(NOT IRL_BUILD_FORTRAN)
  set_target_properties(irl_c irl_fortran 
			PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

# Directory for test files
# Will automatically download GoogleTest and compile it to
# use for the testing framework.
if(${BUILD_TESTING})
  include(CTest)
  enable_testing()
  add_subdirectory(./tests)
endif()

# Directory with source files.
add_subdirectory("${PROJECT_SOURCE_DIR}/irl")

# Directory with examples
set(IRL_EXAMPLES_OUTPUT_DIR "${CMAKE_BINARY_DIR}/compiled_examples")
add_subdirectory("${PROJECT_SOURCE_DIR}/examples")

# Directory with development tools.
#add_subdirectory("${PROJECT_SOURCE_DIR}/tools")

# Export targets to be used by other programs
install(TARGETS irl EXPORT "irl_exp" ARCHIVE DESTINATION "lib")
install(DIRECTORY "${CMAKE_SOURCE_DIR}/irl" # source directory
        DESTINATION "include" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h"# select header files
)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/irl" # source directory
        DESTINATION "include" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.tpp"# select header files
)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/external/abseil-cpp" # source directory
        DESTINATION "include/external" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.h"# select header files
)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/external/abseil-cpp" # source directory
        DESTINATION "include/external" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.inc"# select header files
)

if(IRL_BUILD_FORTRAN)
  install(TARGETS irl_c EXPORT "irl_c_exp" ARCHIVE DESTINATION "lib")
  install(TARGETS irl_fortran EXPORT "irl_fortran_exp" ARCHIVE DESTINATION "lib")
  install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/" DESTINATION "include/irl_fortran")
endif()
