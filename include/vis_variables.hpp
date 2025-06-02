#ifndef HEXED_VIS_VARIABLES_HPP_
#define HEXED_VIS_VARIABLES_HPP_

#include "Namespace.hpp"
#include "Element.hpp"
#include "connection.hpp"

//! \brief functions that assign visualization data to HIL variables
namespace hexed::vis_variables {

/*! \details Assigns the follwing variables:
 * - `is_extruded`: 1 if element is extruded, else 0
 * - `n_dim`: number of dimensions
 * - `is_def = elem.get_is_deformed()`
 * - `ref_level = elem.refinement_level()`
 * - `aniso_ref_level = elem.aniso_ref_level()`
 * - `mask = elem.mask()`
 * - `nom_sz` = elem.nominal_size()`
 * - `uncertainty` = elem.uncertainty`
 * - `sharp`: `true` iff at least one of the element's vertices has been snapped to a sharp edge or point.
 * - `center0`, `center1`, `center2`: center of mass of vertices (not necessarily of the element itself)
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
 */
void surface(Namespace&, Boundary_connection&);

}
#endif
