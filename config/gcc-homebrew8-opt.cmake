set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type")

set(CMAKE_CXX_COMPILER "g++-8" CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER "gfortran-8" CACHE STRING "Fortran Compiler")


set(IRL_CXX_FLAGS "-O3 -DNDEBUG -DNDEBUG_PERF"
    CACHE STRING "C++ compile flags")

set(IRL_Fortran_FLAGS "-O3 -DNDEBUG -DNDEBUG_PERF -ffree-line-length-none"
    CACHE STRING "Fortran compile flags")