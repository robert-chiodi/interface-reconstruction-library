C / Fortran IRL Interface
================

IRL, while written in C++14, contains a C and Fortran2008 interface. Documenting this interface is currently ongoing, but some examples do exist in the `examples/fortran` directory.
It should also be noted, while not thoroughly documented, the Fortran interface has already been used in several large scale CFD flow solvers and is capable of delivering the majority of IRL's functionality.

To use the Fortran IRL interface with Abseil enabled in the build, you must link your application to `libirl_fortran.a`, `libirl_c.a`, `libirl.a`, and `libabsl_all.a`, in that order. The IRL libraries will be in `$install_dir/lib`, and `libabsl_all.a` will be in `$install_dir/absl/lib`. The module files will also need to be included during compilation of your application, and these are placed in `$install_dir/include/irl_fortran`. 
