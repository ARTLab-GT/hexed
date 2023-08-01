<img src="../assets/header.png" alt="header" height="400"/>
Above: Temperature (K) contours on a Mach 10 starfish (inviscid, 2D, CPG).

# hexed
Welcome to Hexed, a Discontinuous Galerkin CFD solver with automated, unstructured quad/hex meshing!

__hex__ (_noun_)
1. (technical) abbreviation for "hexahedron", a six-sided solid (e.g., a cube)
2. curse; jinx; evil magic spell

__hex__  (_verb_)
1. to cast a hex

Anyone who has spent a substantial amount of time working with computational fluid dynamics (CFD)
can attest that it is decidedly arcane and sometimes tends toward evil (in the [programming sense](https://isocpp.org/wiki/faq/big-picture#defn-evil) of the word).
As its name suggests, Hexed scarcely presumes to change that.
What it *does* aim to provide is a more accurate and automated solver, enabling you to practice your witchcraft on a previously unattainable scale.
Specifically, it addresses the following problems with conventional CFD:
- Meshing usually takes a lot of the user's time.
  Although automated meshing algorithms exist, they tend to be fragile when applied to complex geometry.
  This is especially true when attempting to generate quad/hex meshes, on which solvers tend to perform better.
- Solvers are generally slow.
  Most solvers work well for relatively simpler problems but scale poorly (at least, compared to other branches of computational science)
  when the flow involves complex features or high accuracy is required.
- The results are highly dependent on user inputs,
  as there are usually not sufficient computational resources to obtain a fully grid- and time-converged solution.
  User inexperience or insufficient effort can lead to results that look plausible but in fact are completely incorrect.
  For this reason, many engineers distrust CFD.

Hexed (or perhaps "Vexed", as a certain insightful friend suggested it might be more aptly named)
is an all-hex meshing algorithm with a high-order DG solver.
The meshing algorithm has a straightforward geometric interpretation, making it robust to challenging geometries (e.g. small gaps or non-watertight surfaces),
and you can rest assured that whatever output you get results from simple math
and not you, or some sophisticated AI-adjacent algorithm, doing something quirky.
Furthermore, the low-dissipation DG scheme generally gives highly accurate results (when it doesn't blow up),
and when it doesn't, it fails in the form of obviously spurious oscillations.
Thus when you manage to get an answer you can generally be confident in its accuracy.

Hexed is written in C++ and is currently capable of solving the compressible, laminar Navier-Stokes equations for a calorically-perfect gas.
It uses Cartesian meshes with hanging-node refinement, modified to include body-fitted cells at the geometry surface.
In the relatively near future, I hope to include adaptive mesh refinement, moving geometry, and RANS turbulence modeling.
It may not be the most robust or user-friendly, but I hope it makes your projects, if not easier, at least less tedious.

Double, `double`, toil and trouble; CPU burn and contours bubble...

## Contact
\anchor me
This code is developed and maintained by me, Micaiah Smith-Pierce.
I am a PhD student of Professor Stephen Ruffin at the Aerothermodynamics Research and Technology Laboratory (ARTLab),
of the Daniel Guggenheim School of Aerospace Engineering at Georgia Tech.
If you have any questions or feedback, you can reach me at micaiah@gatech.edu or on [Github](https://github.gatech.edu/mcsp3).
