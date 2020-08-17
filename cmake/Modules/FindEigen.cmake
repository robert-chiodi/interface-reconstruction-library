# This module finds the Eigen directory to set an include path to it.
# EIGEN_FOUND is set to True if it is succesfull found.
#
#  EIGEN_INCLUDE_DIR - required include directory
#
# Set the variable EIGEN_PATH to provide a hint to the module for
# where to find the library and header file.  This is searched before the
# standard system locations.
#

find_path(EIGEN_INCLUDE_DIR Eigen
          HINTS "${EIGEN_DIR}")

# Set EXODUS_FOUND to TRUE of the REQUIRED variables have been set.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen
                                  REQUIRED_VARS EIGEN_INCLUDE_DIR)

