# Notation
This document specifies notational conventions used in CartDG. If the meaning
of a number or variable is unclear or ambiguous, look here for guidance. Also, please
try to maintain these conventions when contributing to the code.

## Units and physical quantities
* Unless otherwise specified, all dimensional numbers are in the applicable SI units
  (i.e., derived from m, kg, s).
* Unless otherwise specified, all angles are radian.
* If an extensive thermodynamic quantity (e.g., mass, volume, energy, as opposed to an
  intensive quantity like temperature or pressure) is specified without a quantity, it
  means per unit volume. For example, "energy" specified without other context will be
  in units of J/m^3.
  * As a result, "mass" and "density" are synonymous. In fact, prefer "mass" so that
    all the conserved variables (momentum, mass, energy) are consistently extensive.
* Stagnation quantities will be referred to as such. The word "total" shall not be
  used to indicate stagnation quantities.
* The flow state vector (often simply referred to as `state`) consists of the
  conserved variables.
  It has components [momentum, mass, total energy]. E.g., in 3D that is
  [momentum0 (x-direction), momentum1 (y-direction), momentum2 (z-direction), mass, energy].
  "Total energy" means thermodynamic internal energy + kinetic energy.
  * This is different from most codes, which use the order (mass, momentum, energy). The
    order for this code was chosen because it groups the vector and scalar quantities
    together and because the i'th component of the state vector is equal to the i'th
    component of momentum. The latter makes the code for certain formulae cleaner.

## Storage order
* Multidimensional arrays are row-major by default. However, be careful, because
  Eigen matrices are column-major by default, so working with them often involves
  transposing the data.
* When storing state data in an array, the order of the indices is usually the
  following: [Element, Runge-Kutta stage, state component, x-index, (y-index, (z-index))].
* Vertices are ordered as a 2(x2(x2)) array.
* When storing face data as an array, the order of indices is usually:
  [dimension-index, positivity (negative side comes before positive), (x-index), (y-index), (z-index)].
  Of course, for 3D arrays at most two of the x, y, and z indices are included, and for 2D at most two.

## Abbreviations
In some abbreviations that are used very frequently, clarity is traded for brevity.
Thus, please be familiar with the following definitions:
* `i_`: "index of ..."
  * If this is an index in a nested loop, the inside loops may use `j_` and `k_`.
* `n_`: "number of ..."
* `qpoint`: "quadrature point"
* `heat_rat`: "ratio of heat capacities". Often denoted "gamma" in thermodynamics.
* `dim`: dimension
  * `i_dim` specifically serves to identify which of the coordinate axes is referenced
    in some operation.
* `def`: "deformed"
* `con`: "connection"
* `ener`: "total energy"
* `_sq`: "... squared"

## Terminology
CartDG uses some terms that are non-standard, or at least not completely universal.
They are defined as follows:
* "kernel": A performance-critical function which utilizes template metaprogramming
  and relies on a wrapper function to select the appropriate template based on runtime
  parameters.
* "local" kernel: performs updates based on information stored in a single element.
* "neighbor" kernel: calculates adjustments based on numerical flux at the interface
  between two elements
* The "row size" of a basis is its degree + 1 (i.e., the number of rows of quadrature
  points in each dimension, which equals the number of coefficients in the 1D case).
  It is usually more convenient to talk about the row size than the degree. The row
  size is not to be confused with the number of quadrature points, which depends on
  the dimensionality. E.g., a 4th degree basis in 3 dimensions has row size 5 and 125
  quadrature points.
* "node": An interpolation node for a basis of Lagrange (interpolating) polynomials.
  These tend to coincide with the quadrature points.
* "vertex": A single point shared by a collection of elements. Not to be confused
  with "node".
