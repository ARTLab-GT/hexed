# hexed

__hex__ (_noun_)
1. (technical) abbreviation for "hexahedron", a six-sided solid (e.g., a cube)
2. curse; jinx; evil magic spell

__hex__  (_verb_)
1. to cast a hex

[full documentation page](https://github.gatech.edu/pages/ARTLab/hexed/)

Hexed is a Discontinuous Galerkin CFD solver with automated, unstructured quad/hex meshing.
It is designed to address the following problems with conventional CFD:
- Meshing usually takes a lot of the user's time.
  Although automated meshing algorithms exist, they tend to be fragile when applied to complex geometry.
  This is especially true when attempting to generate quad/hex meshes, on which solvers tend to perform better.
- Solvers are generally slow.
  Most solvers work well for relatively simpler problems but scale poorly (at least, compared to other branches of computational science)
  when the flow involves complex features or high accuracy is required.
- The results are highly dependent on user inputs,
  as there are usually not sufficient computational resources to obtain a fully grid- and time-converged solution.
  User inexperience can lead to results that look plausible but in fact are completely incorrect.
  For this reason, many engineers distrust CFD.

Hexed's all-hex meshing algorithm has a straightforward geometric interpretation,
making it robust to challenging geometries (e.g. small gaps or non-watertight surfaces),
and you can rest assured that whatever output you get results from simple math and not you,
or some sophisticated AI-adjacent algorithm, making a questionable decision.
Furthermore, the low-dissipation DG scheme generally gives highly accurate results (when it doesn't crash),
and when it doesn't, it fails in the form of obviously spurious oscillations.
Thus when you manage to get an answer you can generally be confident in its accuracy.

Hexed is written in C++ and is currently capable of solving the compressible, laminar Navier-Stokes equations for a calorically-perfect gas.
It uses Cartesian meshes with hanging-node refinement, modified to include body-fitted cells at the geometry surface.
In future versions, I hope to include adaptive mesh refinement, moving geometry, and RANS turbulence modeling.
