set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build Type")

set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "Fortran Compiler")


set(IRL_CXX_FLAGS "-g -O0 -debug all -ftrapuv -traceback -pedantic -Wall -Wextra"
    CACHE STRING "C++ compile flags")

set(IRL_Fortran_FLAGS "-g -traceback -check all -check bounds -check uninit -ftrapuv"
    CACHE STRING "Fortran compile flags")
