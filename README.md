<img src="../assets/header.png" alt="header" height="400"/>
Above: Temperature (K) contours on a Mach 10 starfish (inviscid, 2D, CPG).

# hexed
Discontinuous Galerkin engine for CFD with automated, unstructured quad/hex meshing.

**hex** *(noun)*
-# (technical) abbreviation for "hexahedron", a six-sided solid (e.g., a cube)
-# curse; jinx; evil magic spell

**hex** *(verb)*
-# to cast a hex

This document provides an overview of Hexed. More detailed documentation can be found in [`doc`](doc/). In particular, for
installation instructions, see [`install.md`](doc/install.md).

## Overview
Anyone who has spent a substantial amount of time working with computational fluid dynamics (CFD)
can attest that it is decidedly arcane and sometimes tends toward evil.
As its name suggests, Hexed scarcely presumes to change that.
What it *does* aim to provide is a faster, more automated solver, enabling you to practice your witchcraft on a previously unattainable scale.
Specifically, Hexed (or perhaps "Vexed", as a certain insightful friend suggested it might be more aptly named)
is a C++ library which can solve the compressible Navier-Stokes equations of aerodynamics on unstructured quad/hex meshes.
It can handle Cartesian meshes with hanging-node refinement, and it also provides a mechanism to automatically generate a body-fitted mesh
from a Cartesian starting point.
The high-order discontinuous Galerkin (DG) scheme is designed to achieve whatever level of accuracy is required with a much coarser mesh,
reducing (hopefully) the overall computational cost.
Hexed is used as a backend by [NASCART-GT](https://github.gatech.edu/ARTLab/NASCART-GT)
to combine the speed and accuracy of the DG method with the automation and versatility of its
adaptive Cartesian mesh.
However, it is meant to be encapsulated well enough that it could survive a complete refactor of NASCART-GT,
or even be used with other frontends.
Hexed is still very much a work in progress, but hopefully it can make your CFD, if not easier, at least less tedious.

To summarize, Hexed is:
* fast
* high-order accurate
* automatable
* unstructured

Hexed is not:
* a stand-alone application
* for arbitrary PDEs
* for tri/tet grids

## Features
Currently implemented features:
* Unstructured, body-fitted quad/hex mesh.
* Explicit time integration.
* Visualization with Tecplot.
* Integration with NASCART-GT.
* Hanging-node refinement.
* Shock capturing with an original artificial viscosity scheme.
* Viscous flows with anisotropic refinement.

Planned or in-progress features (roughly in order of planned implementation):
* Grid adaption.
