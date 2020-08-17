set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type")

set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "Fortran Compiler")


set(IRL_CXX_FLAGS "-O3 -DNDEBUG -DNDEBUG_PERF -fPIC"
    CACHE STRING "C++ compile flags")

set(IRL_Fortran_FLAGS "-O3 -DNDEBIG -DNDEBUG_PERF -fPIC"
    CACHE STRING "Fortran compile flags")
