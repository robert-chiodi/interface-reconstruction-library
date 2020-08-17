set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build Type")

set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "Fortran Compiler")


set(IRL_CXX_FLAGS "-g -O2 -DNDEBUG -DNDEBUG_PERF -fPIC"
    CACHE STRING "C++ compile flags")

set(IRL_Fortran_FLAGS "-g -O2 -DNDEBUG -DNDEBUG_PERF -fPIC"
    CACHE STRING "Fortran compile flags")
