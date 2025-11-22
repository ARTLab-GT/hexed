#ifndef HEXED_VIS_VARIABLES_HPP_
#define HEXED_VIS_VARIABLES_HPP_

#include "Namespace.hpp"
#include "Element.hpp"
#include "Boundary_connection.hpp"

//! \brief functions that assign visualization data to HIL variables
namespace hexed::vis_variables {

/*! \details Assigns the follwing variables:
 * - `is_extruded`: 1 if element is extruded, else 0
 * - `n_dim`: number of dimensions
 * - `is_def = Element::get_is_deformed()`
 * - `ref_level = Element::refinement_level()`
 * - `aniso_ref_level = Element::aniso_ref_level()` __deprecated__
 * - `aniso_ref_level0`, `aniso_ref_level1`, `aniso_ref_level2`: The `Tree::anisotropic_refinement_level`
 *   of this element's `Tree`.
 * \warning For now, `aniso_ref_level` and `aniso_ref_level0`, ... do completely different things.
 * - `mask = Element::mask()`
 * - `nom_sz = Element::nominal_size()` __deprecated__
 * - `nominal_size = Element::nominal_size()`
 * - `wall_distance = Element::wall_distance()`
 * - `wall_dimension = Element::wall_dimension()`
 * - `has_wall = Element::has_wall()`
 * - `uncertainty = Element::uncertainty`
 * - `snapping_problem`: Whether there was any problem snapping this element to geometry features.
 * - `is_deformed`: True if this element is deformed, false if it is Cartesian.
 * - `sharp`: `true` iff at least one of the element's vertices has been snapped to a sharp edge or point.
 * - `center0`, `center1`, `center2`: center of mass of vertices (not necessarily of the element itself)
 * - `nominal_shape0`, `nominal_shape1`, `nominal_shape2`: `Element::nominal_shape(i_dim)` for `i_dim` in [0, 3).
 *   Trailing dimensions set to 0.
 * - `rms_residual`: RMS (within the element, in reference space) of the residual of all flow variables in the element.
 *   Residuals of all variables are summed.
 * - `spectral_uncertainty0`, `spectral_uncertainty1`, `spectral_uncertainty2`:
 *   An anisotropic measure of the spectral convergence of the flow variables in this element.
 *   For each variable and each dimension,
 *   the magnitude of the highest-order polynomial mode in each dimension normalized by the difference between
 *   the maximum and the minimum of that variable over the entire domain is computed.
 *   The spectral uncertainty is then set to the maximum of this metric for all mean flow variables
 *   (but not turbulence variables).
 * - `flux_uncertainty`: For wall elements, the norm of the jump in the momentum flux between this element and
 *   the next farthest element from the wall normalized by the norm of the momentum flux on that same face.
 *   For non-wall elements, 0.
 *   \todo Residual contributions from each variable need to be normalized.
 *   Right now it might as well be the energy residual.
 * - `diffusion_limited`: 1 if the time step constraint imposed by diffusion is stricter
 *   than the time step constraint imposed by convection at any of the quadrature points, and 0 otherwise.
 * - `source_limited`: 1 if the time step constraint imposed by source_terms is stricter
 *   than the time step constraint imposed by convection and diffusion combined
 *   at any of the quadrature points, and 0 otherwise.
 *
 * \deprecated The following subset of the assigned variables are deprecated and will be removed in a future version:
 * \deprecated
 * `nom_sz`: Use the synonymous `nominal_size` variable instead.
 * \deprecated
 * `aniso_ref_level`: This referred to anisotropy produced by \ref split_layers, which itself is deprecated.
 * To inspect the anisotropy produced by anisotropic refinement, use `aniso_ref_level0`, `aniso_ref_level1`,
 * `aniso_ref_level2`.
 */
void element(Namespace&, Element& elem);

/*! \brief Assigns the variables `pos0`, `pos1`, `pos2` and `jacobian_det`.
 * \details to the position of the `i_qpoint`th quadrature point and the Jacobian determinant.
 * Trailing dimensions are set to 0.
 */
void position(Namespace&, Element&, const Basis&);

/*! \details Assigns the follwing variables:
 * - `momentum0`, `momentum1`, `momentum2` : momentum per volume
 * - `mass`: mass per volume (aka density)
 * - `energy`: total energy per volume
 * - `tss`: time step scale
 * - `art_visc`: artificial viscosity coefficient
 * - `residual_xxx` for `xxx` in {`momentum0`, ..., `energy`}: residual of each of the conserved state variables
 *
 * If turbulent, also assigns the following variables:
 * - `turbulent_kinetic_energy`: the \f$ \rho k \f$ in two-equation turbulence models
 * - `turbulent_dissipation_bassi`: \f$ \rho \tilde{\omega} = \rho \ln \omega \f$ in the \f$ k\text{-}\omega \f$ model.
 *   Note that this technically violates the rules of dimensional analysis
 *   by taking the log of a dimensional quantity,
 *   but numerically this will not cause a problem.
 *   A change of units will simply manifest as a constant offset on \f$ \tilde{\omega} \f$.
 *
 * \attention If `Solver::update` or `Solver::update_art_visc_smoothness` have been called
 * since the last call to `Solver::compute_residuals` then the residual variables will be incorrect.
 */
void state(Namespace&, Element&);

//! \brief Assigns everything in `element()`, `position()`, and `state()`.
void field(Namespace&, Element&, const Basis&);

/*! \details Assigns the follwing variables:
 * - `pos0`, `pos`, `pos2`: position
 * - `normal0`, `normal1`, `normal2`: unit surface normal (out of surface, into domain)
 * - `momentum0`, `momentum1`, `momentum2` : momentum per volume
 * - `mass`: mass per volume (aka density)
 * - `energy`: total energy per volume
 * - `visc_stress0`, `visc_stress1`, `visc_stress2` : viscous stress at surface
 * - `mass_flux`: diffusive mass flux through surface
 * - `heat_flux`: surface heat flux
 * - `roughness_height`: equivalent sand roughness height (dimensional) used for \f$ \omega \f$ boundary condition
 * - `wall_spacing`: wall distance of first element's farthest vertex from the wall
 */
void surface(Namespace&, Boundary_connection&);

}
#endif
