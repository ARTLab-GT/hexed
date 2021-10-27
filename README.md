# CartDG
Provides a discontinuous Galerkin scheme as a library for Cartesian grid CFD solvers.

## Overview
CartDG implements the core numerical algorithms for a discontinuous Galerkin scheme operating on the equations of fluid flow
defined on a Cartesian grid. It takes advantage of mathematical simplifications resulting from the Cartesian,
isotropic nature of the majority of the grid to produce concise and efficient code. It is designed to be integrated with NASCART-GT as
an accuracy boost, but it is meant to be encapsulated well enough that it can survive major changes to NASCART-GT or even be
used with other similar codes.

CartDG is:
* An implementation of the discontinuous Galerkin scheme for CFD solvers.
* Fast.
* Simple (relatively).

CartDG is not:
* A stand-alone CFD code.
* A general framework for CFD on arbitrary grids.
* A solver for arbitrary PDEs.

This document describes the purpose and capabilities of CartDG. For installation instructions, see [`INSTALL.md`](INSTALL.md).
For contributing guidelines, see [`CONTRIBUTING.md`](CONTRIBUTING.md).

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
The speed of the kernels is measured by the script [`benchmark.py`](script/benchmark.py). The speed of the code as a whole has shown good agreement with the measurements made by `benchmark.py`.

## Terminology and conventions
I use the following nonstandard terms in the code:
* "local" kernel - performs updates based on information stored in a single element.
* "neighbor" kernel - the kernel that calculates adjustments based on numerical flux at the interface between two elements.
* The "row size" of a basis is its degree + 1 (i.e. the number of rows of quadrature points in each dimension, which equals
  the number of coefficients in the 1D case). It is usually more convenient to talk about the row size
  than the degree. The row size is not to be confused with the number of quadrature points, which depends on the
  dimensionality. E.g. a 4th degree basis in 3 dimensions has row size 5 and 125 quadrature points.
 
## Dependencies
Eigen must be available in your include path. Tecplot must be installed with path environment variables configured accordingly.
Catch2 must be available in a location that CMake can find.
Python3 must also be available, along with the libraries NumPy, SymPy, and MatPlotLib. See [installation instructions](INSTALL.md)
for guidance on obtaining these.
 
## Features
Currently implemented features:
* Standard DG method.
* Gauss-Lobatto and Gauss-Legendre nodal bases.
* Body-fitted quasi-Cartesian quad mesh.
* Cartesian hex mesh.
* 3-stage Runge-Kutta explicit time integration.
* CFL-based time step selection.
* Visualization with Tecplot.
* Integration with NASCART-GT.
 
Planned or in-progress features (roughly in order of planned implementation):
* Shock capturing (in progress).
* Body-fitted quasi-Cartesian hex mesh.
* Isotropic hanging-node refinement.
* Viscous flows with anisotropic refinement.
* Grid adaptation.
