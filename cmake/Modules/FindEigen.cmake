#[=======================================================================[.rst:
FindEigen
-------

Finds the Eigen library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``Eigen::Eigen``
The Eigen library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``Eigen_FOUND``
True if the system has the Eigen library.
``Eigen_VERSION``
The version of the Eigen library which was found.
``Eigen_INCLUDE_DIRS``
Include directories needed to use Eigen.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Eigen_INCLUDE_DIR``
The directory containing ``Eigen.h``.Eigen

#]=======================================================================]

# Try to grab information from PkgConfig
find_package(PkgConfig)
pkg_check_modules(PC_Eigen QUIET Eigen)

if(NOT EIGEN_DIR)
  find_path(EIGEN_INCLUDE_DIR
  NAMES Eigen/Dense
  PATHS ${PC_EIGEN_INCLUDE_DIRS}
  )
else()
  find_path(EIGEN_INCLUDE_DIR
  NAMES Eigen/Dense
  PATHS ${EIGEN_DIR}
  NO_DEFAULT_PATH
  )
endif()

set(Eigen_VERSION ${PC_Eigen_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen
FOUND_VAR EIGEN_FOUND
REQUIRED_VARS
EIGEN_INCLUDE_DIR
VERSION_VAR EIGEN_VERSION
)

if(NOT EIGEN_FOUND)
message("Eigen not found. Define EIGEN_DIR with path to Eigen.")
return()
endif()

if(EIGEN_FOUND)
set(EIGEN_INCLUDE_DIRS ${Eigen_INCLUDE_DIR})
set(EIGEN_DEFINITIONS ${PC_Eigen_CFLAGS_OTHER})
endif()

if(EIGEN_FOUND AND NOT TARGET Eigen::Eigen)
add_library(Eigen::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen::Eigen PROPERTIES
INTERFACE_COMPILE_OPTIONS "${PC_Eigen_CFLAGS_OTHER}"
INTERFACE_INCLUDE_DIRECTORY "${EIGEN_INCLUDE_DIR}"
INTERFACE_INCLUDE_DIRECTORIES "${EIGEN_INCLUDE_DIR}"
)
endif()

mark_as_advanced(
EIGEN_INCLUDE_DIR
)

