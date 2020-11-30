Available Examples
===================================

For available examples, see the `examples` directory. The examples consist of some basic C++ and Fortran capabilties.
Additionally, in the `examples/simple_advector` and `examples/advector`, complete geometric VOF advection algorithms are implemented. 

In order to run the examples, IRL must first be configured using CMake. See the installation documentation for directions on how to do that. 
Afterwards, the C++ and Fortran examples can be compiled using the Make targets `examples_cpp` and `examples_fortran`, respectively.
The advector examples can similarly be built with the Make targets `examples_simple_advector` or `examples_advector`.
Upoon being built, these examples will be placed in the build directory in the `compiled_examples` directory.
