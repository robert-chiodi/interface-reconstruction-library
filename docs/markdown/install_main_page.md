# How to Install and Link to IRL

IRL uses CMake (version >= 3.11) to handle its configuration and compilation. By default, the IRL C++, C, and Fortran interfaces will be built and placed in the `install` directory within a `Release`, `Debug`, or `RelWithDebInfo` (depending on the type of build requested). After installing the dependencies (detailed below), the simplest IRL CMake build command can be executed as
```
mkdir build &&
cd build &&
cmake -C ../config/gcc-opt.cmake -D EIGEN_DIR=/path/to/Eigen .. &&
make install
```
where it has been assumed that GCC compilers will be used and your own path to Eigen has been substituted for `/path/to/Eigen`. 

Changing the compiler used or the compiler options passed to IRL is straight-forward, and just requires creating and modifying a new *.cmake file in the `config` directory. In `config` there exists many examples, one of which will most likely already work for you.

Several additional options can be passed to the CMake configuration to tailor IRL to your needs. They are described below.

* `-D EIGEN_DIR=/path/to/Eigen` ----- Provides path to Eigen directory.
* `-D BUILD_TESTING=ON` ----- Turns on building of unit tests, which are not built by default.  If building tests, GoogleTest will be downloaded and compiled unless the GOOGLETEST_DIR CMake flag is given. 
* `-D GOOGLETEST_DIR=/path/to/google_test` ----- Provides path to [Google Test](https://github.com/google/googletest) that will be used if when building tests (if turned on).
* `-D USE_ABSL=OFF` ----- IRL ships with a version of [Google Abseil](https://github.com/abseil/abseil-cpp) for use of its `inlined_vector` and `flat_hash_map` classes. By turning it off, IRL will revert back to its own implementation.


## External Dependencies

### Eigen

The header-only matrix library [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is used inside IRL. Before configuring IRL, Eigen should be downloaded from their website into a directory of your choice. During CMake configuration of IRL, this path is then supplied to IRL via `-D EIGEN_DIR=/path/to/Eigen`.

### Google Test

If building IRL's unit tests, the [Google Test](https://github.com/google/googletest) framework is needed. By default, Google Test will be downloaded during CMake configuration and compiled alongside the building of IRL's unit tests. If Google Test is already installed locally, it can be used instead via the CMake variable `-D GOOGLETEST_DIR=/path/to/google_test`. This is also required for machines that cannot readily access GitHub to download Google Test directly.
