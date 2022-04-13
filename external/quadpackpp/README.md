# Quadpack++ adaptive quadrature

_The software is no longer active development. The author periodically
receives emails asking about its availability, and has decided to make GitHub
its permanent home. Reasonable effort will be made to address open issues._

## What it does

quadpack++ implements adaptive quadrature routines from the
[GNU Scientific Library](http://www.gnu.org/software/gsl/)
for arbitrary precision floating-point arithmetic by using C++ templates.
The core routines are implemented in header files, requiring nothing to be
separately built and linked against.

## How to use it

The program [examples/logarithmic.cpp](examples/logarithmic.cpp) illustrates
how to use the library to integrate the function _f(x) = x^α log(1/x)_ on the
unit interval. (This function is part of a standard test battery for numerical
quadrature methods. It is integrable for all values α > 0, but is challenging
for non-adaptive routines near the left endpoint for parameter values that
allow the logarithmic singularity to dominate.)

For floating-point types with precision beyond native (eg, `long double`),
[MPFR C++](http://www.holoborodko.com/pavel/mpfr/) has been successfully used
and appears to receive continuous maintenance support.

## Design principle

The GSL routines (like QAG, QAGS, etc) on which this library is modeled
require separate "workspace" to be allocated that is passed as an argument
to the static methods. Here, [Workspace](include/workspace.hpp) is a class
that instantiated with the scratch-space size and base "degree" arguments,
and all quadrature methods are exposed as class members.

Like the GSL routines, quadrature methods take as arguments either absolute
or relative error tolerances, and they return both the result and a numerical
error estimate.

### Gauss-Kronrod quadrature

The Gauss-Kronrad quadrature algorithm makes the robust error estimation
and control possible. This algorithm relies on a set of precomputed abscissae
(in the unit interval) and corresponding weights for a given degree _2m+1_.
In the original QUADPACK/GSL routines these were tabulated as constants in
the source code only for cases m=7, 10, 15, 20, 25, and 30. Here the
[GaussKronrod](include/gauss-kronrod.hpp) class dynamically computes the
abscissae and weights, creating a GK-rule for any valid _m ≥ 1_ and any
floating-point precision templated by the `Real` type.

Error control in the adaptive routines relies on "machine" parameters
associated to the floating point type, specifically
[machine epsilon](https://en.wikipedia.org/wiki/Machine_epsilon),
overflow limit, and underflow limit. These are dynamically computed in
the [Machar](include/machar.hpp) class, which is based on the
[MACHAR](https://people.sc.fsu.edu/~jburkardt/f77_src/machar/machar.f)
routine original published in ACM Transactions of Mathematical Software.
