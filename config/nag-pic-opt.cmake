set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type")

set(CMAKE_CXX_COMPILER "g++" CACHE STRING "C++ Compiler")
set(CMAKE_Fortran_COMPILER "nagfor" CACHE STRING "Fortran Compiler")


set(IRL_CXX_FLAGS "-O3 -DNDEBUG -DNDEBUG_PERF -fPIC"
    CACHE STRING "C++ compile flags")

set(IRL_Fortran_FLAGS "-O3 -DNDEBUG -DNDEBUG_PERF -PIC"
    CACHE STRING "Fortran compile flags")

