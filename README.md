# CartDG
Discontinuous Galerkin numerics for CFD on a Cartesian grid.

## Overview
CartDG implements the core numerical algorithms for a discontinuous Galerkin scheme operating on the equations of fluid flow
defined on a Cartesian grid. It takes advantage of mathematical simplifications resulting from the Cartesian,
isotropic nature of the grid to produce concise and efficient code. It is designed to be integrated with NASCART-GT as
an accuracy boost, but it is meant to be encapsulated well enough that it can survive major changes to NASCART-GT or even be
used with other similar codes.

What CartDG is:
* Fast.
* Simple.
* Capable of helping NASCART-GT do CFD with a state-of-the art numerical scheme.

What CartDG is not:
* A stand-alone CFD code.
* A general framework for CFD on arbitrary grids.
* Capable of solving equations other than those that govern the flow of fluids.

## Motivation
The motivation for this project came from my experience with NASCART-GT. With it's highly automated adaptive Cartesian grid
strategy, NASCART-GT represents the holy grail of modern CFD, in that the user merely provides the geometry and a few
basic parameters and the code does the rest. However, it still takes a long time to solve even relatively basic problems
(which is mostly a feature of mainstream CFD methods rather than of NASCART-GT in particular). For that reason I am attempting
to use a high-order-accurate discontinuous Galerkin method to achieve the same accuracy with many fewer cells and less overall
time.

## Some notes on the implementation
Most of the performance critical code is placed in what I call "kernels". These are highly optimized functions that perform
the basic operations of the DG method. They take things like the number of dimensions and the polynomial degree of the basis
as template parameters, allowing Eigen matrices and local arrays to be allocated on the stack and allowing the compiler to
unroll loops effectively. The ability to specify the degree of the basis and the dimensionality at runtime is recovered by
a very ugly file that assigns an array of function pointers to kernels with different template arguments (this file is
generated automatically by a python script). The kernels also make minimal use of abstraction. This means that the same
code is sometimes written multiple times, but it gives the compiler maximal freedom to optimize.

All of these things are
disasterous from a software engineering perspective, but have shown significant performance benefits. For the kernels
I am willing to make this trade, but for the rest of the code I have priortized readability and modularity over performance.
The speed of the kernels is measured by the script `benchmark.py`. The speed of the code as a whole when solving the
isentropic vortex problem has shown good agreement with the measurements made by `benchmark.py`.

## Terminology and conventions
I use the following nonstandard terms in the code:
* "local" kernel - the kernel that calculates time derivatives based only on information within each element.
* "neighbor" kernel - the kernel that calculates adjustments based on numerical flux at the intersection of (the closures
   of) two elements.
* At some point in the future I may begin using the term "parent" and "child" kernel to refer to kernels that restrict and
  prolong the solution onto coarser or finer grids in a multigrid sense.
* The "rank" of a basis is it's degree + 1 (i.e. the number of coefficients in it's 1-D polynomial). This is generally a
  more convenient term to use than the degree, and it is not to be confused with the number of quadrature points, which
  depends on the dimensionality. E.g. a 4th degree basis in 2 dimensions is of rank 5 and has 25 quadrature points. This
  unfortunately may conflict with the use of the term "rank" in parallel programming, but I have yet to invent a better
  term.
 
## Compiling
Perform the following steps to build the library `libcartdg.so` which contains most of the functionality.
* Generate translation files that associate function pointers with instances of kernel templates and define
  high-precision quadrature rules.
  * `cd script`
  * `python3 auto_generate.py`
  * `cd ..`
* Edit the
  Makefile to select the appropriate compiler flags. Use the debug flags (the ones with `-fsanitize`) if you want to run the
  tests. Use the high-performance flags (the ones with `-O3` and `-march=native`) if you wish to run the demo.
* `make`

To build and run the tests:
* `cd test`
* `make`
* `./test`

All tests should pass, and there should
also be some Tecplot files.

To build and run the demo:
* `cd demo`
* `make`
* `./demo`

The executable will simulate the isentropic vortex test problem and write the solution to `.szplt` Tecplot binary files.
 
This is obviously far from an ideal build system, but it works for now until I learn how to properly set up CMake.

## Dependencies
The following must be available in your include and library paths:
* Eigen
* Tecio
 
## Features
Currently implemented features:
* Standard DG method.
* Gauss-Lobatto nodal basis.
* Isotropic, Cartesian cells connected in an arbitrary graph.
* 3-stage Runge-Kutta explicit time integration.
* CFL-based time step selection.
* Generation of structured grids for testing purposes.
* Visualization with Tecplot.
* Periodic boundary conditions.
 
Planned or in-progress features (roughly in order of planned implementation):
* Entropy-stable methods (in progress).
* Integration with NASCART-GT.
* Immersed and freestream boundary conditions.
* Hanging-nodes.
* Shock capturing.
* Implicit and/or multigrid method.

## Contributing
Please do not push directly to `master`. Feel free to create and push new branches. Once you have something that
works, please open a pull request.
