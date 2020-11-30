# How to Install and Link to IRL

IRL uses CMake (version >= 3.9) to handle its configuration and compilation. By default, the IRL C++, C, and Fortran interfaces will be built. The install directory can be modified in typical cmake fashion by defining `-D CMAKE_INSTALL_PREFIX=install_directory`.
After installing the dependencies (detailed below), the simplest IRL CMake build command can be executed as
```
mkdir build &&
cd build &&
cmake .. &&
make install
```
The path to Eigen (a required dependency) can be given through the CMake variable `-D EIGEN_DIR=path/to/Eigen`.  If not supplied, Eigen will be searched for in the CMake default locations, such as `/usr/include`.

To change the compiler or compiler options, change the environment variables `CXX` and `FC` for C++ and Fortran, respectively. The compiler flags can also be changed via the environment variables `CXXFLAGS` and `FFLAGS`. As opposed to settign the environment variables, the above information can also be supplied via the CMake variables `CMAKE_CXX_COMPILER`, `CMAKE_Fortran_COMPILER`, `CMAKE_CXX_FLAGS`, and `CMAKE_Fortran_FLAGS`. Note: If the `CMAKE_BUILD_TYPE` is not `Debug`, the compiler flags `-DNDEBUG` and `-DNDEBUG_PERF` will be appended to `CMAKE_CXX_FLAGS`. The default value for `CMAKE_BUILD_TYPE` is `Release` if one is not supplied.

Several additional options can be passed to the CMake configuration to tailor IRL to your needs. They are described below.

* `-D EIGEN_DIR=/path/to/Eigen` ----- Provides path to Eigen directory.
* `-D BUILD_TESTING=ON` ----- Turns on building of unit tests, which are not built by default.  If building tests, GoogleTest will be downloaded and compiled unless the GOOGLETEST_DIR CMake flag is given. 
* `-D GOOGLETEST_DIR=/path/to/google_test` ----- Provides path to [Google Test](https://github.com/google/googletest) that will be used if when building tests (if turned on).
* `-D USE_ABSL=OFF` ----- IRL ships with a version of [Google Abseil](https://github.com/abseil/abseil-cpp) for use of its `inlined_vector` and `flat_hash_map` classes. By turning it off, IRL will revert back to its own implementation. This is strongly discouraged for performance reasons.


## External Dependencies

### Eigen

The header-only matrix library [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is used inside IRL. Before configuring IRL, Eigen should be downloaded from their website into a directory of your choice. During CMake configuration of IRL, this path is then supplied to IRL via `-D EIGEN_DIR=/path/to/Eigen`.

### Google Test

If building IRL's unit tests, the [Google Test](https://github.com/google/googletest) framework is needed. By default, Google Test will be downloaded during CMake configuration and compiled alongside the building of IRL's unit tests. If Google Test is already installed locally, it can be used instead via the CMake variable `-D GOOGLETEST_DIR=/path/to/google_test`. This is also required for machines that cannot readily access GitHub to download Google Test directly.
