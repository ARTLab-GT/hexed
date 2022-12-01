<img src="../assets/header.png" alt="header" height="400"/>
Above: Temperature (K) contours on a Mach 10 starfish (inviscid, 2D, CPG).

# hexed
Discontinuous Galerkin engine for CFD with automated, unstructured quad/hex meshing.

**hex** *(noun)*
1. (technical) abbreviation for "hexahedron", a six-sided solid (e.g., a cube)
2. curse; jinx; evil magic spell

**hex** *(verb)*
1. to cast a hex

This document provides an overview of Hexed. More detailed documentation can be found in [`doc`](doc/). In particular, for
installation instructions, see [`install.md`](doc/install.md).

## Overview
Anyone who has spent a substantial amount of time working with computational fluid dynamics (CFD)
can attest that it is decidedly arcane and sometimes tends toward evil.
As its name suggests, Hexed scarcely presumes to change that.
What it *does* aim to provide is a faster, more automated solver, enabling you to practice your witchcraft on a previously unattainable scale.
Specifically, Hexed is a C++ library which can solve the compressible Euler equations of aerodynamics on unstructured quad/hex meshes.
It can handle Cartesian meshes with hanging-node refinement, and it also provides a mechanism to automatically generate a body-fitted mesh
from a Cartesian starting point.
The high-order discontinuous Galerkin (DG) scheme is designed to achieve whatever level of accuracy is required with a much coarser mesh,
reducing (hopefully) the overall computational cost.
Hexed is used as a backend by NASCART-GT to combine the speed and accuracy of the DG method with the automation and versatility of its
adaptive Cartesian mesh.
However, it is meant to be encapsulated well enough that it could survive a complete refactor of NASCART-GT,
or even be used with other frontends.
Hexed is still a work in progress and there are many more features to implement, but hopefully it can make your CFD, if not easier, at least less tedious.

To summarize, Hexed is:
* fast
* high-order accurate
* automatable
* unstructured

Hexed is not:
* a stand-alone application
* for arbitrary PDEs
* for tri/tet grids

## Some notes on the implementation
Most of the performance critical code is placed in what I call "kernels".
These are highly optimized functions that perform the basic operations of the DG method.
They take the number of dimensions and the polynomial degree of the basis as template parameters,
providing as many opportunities as possible for compiler optimization (and making it roughly twice as fast).
The ability to specify the degree of the basis and the dimensionality at runtime is recovered by
some truly cursed recursive templates, which gets the end result that the user needs, but takes *forever* to compile
and gives me a headache every time I look at it.
The kernels also often trade abstraction for optimizability.
All of these things are disasterous from a software engineering perspective, but have shown significant performance benefits.
For the kernels I am willing to make this trade, but for the rest of the code I have priortized readability and modularity over performance.
If you end up attempting to understand the kernels, for one reason or another, all I can say is... *sorry*.
 
## Dependencies
The following are required to be installed before compiling Hexed:
- GCC (of course)
- CMake
- Eigen
- Python3
- Numpy
- Scipy
- Sympy

The following are techincally not necessary, but are required by certain optional features:
- Catch2 (unit testing)
- Tecplot (flow visualization, required in order to build with NASCART-GT)
- Otter (flow visualization, not yet publicly available)

See [installation instructions](doc/install.md)
for guidance on obtaining these.
 
## Features
Currently implemented features:
* Gauss-Lobatto and Gauss-Legendre nodal bases.
* Unstructured, body-fitted quad/hex mesh.
* Explicit time integration.
* Visualization with Tecplot.
* Integration with NASCART-GT.
* Hanging-node refinement.
* Shock capturing with an original artificial viscosity scheme.
 
Planned or in-progress features (roughly in order of planned implementation):
* Viscous flows with anisotropic refinement.
* Grid adaption.
