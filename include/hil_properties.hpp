#ifndef HEXED_HIL_PROPERTIES_HPP_
#define HEXED_HIL_PROPERTIES_HPP_

#include "Namespace.hpp"
#include "Element.hpp"
#include "connection.hpp"

//! \brief functions that assign properties of objects to HIL variables
namespace hexed::hil_properties
{

/*! Assigns the follwing variables:
 * - `is_extruded`: 1 if element is extruded, else 0
 * - `n_dim`: number of dimensions
 * - `is_def = elem.get_is_deformed()`
 * - `ref_level = elem.refinement_level()`
 * - `nom_sz` = elem.nominal_size()`
 * - `uncertainty` = elem.uncertainty`
 * - `center0`, `center1`, `center2`: center of mass of vertices (not necessarily of the element itself)
 */
void element(Namespace&, Element& elem);

/*! Assigns the variables `pos0`, `pos1`, `pos2`
 * to the position of the `i_qpoint`th quadrature point.
 * Trailing dimensions are set to 0.
 */
void position(Namespace&, Element&, const Basis&, int i_qpoint);

/*! Assigns the follwing variables:
 * - `momentum0`, `momentum1`, `momentum2` : momentum per volume
 * - `mass`: mass per volume (aka density)
 * - `energy`: total energy per volume
 * - `tss`: time step scale
 * - `art_visc`: artificial viscosity coefficient
 */
void state(Namespace&, Element&, int i_qpoint);

/*! Assigns the follwing variables:
 * - `pos0`, `pos`, `pos2`: position
 * - `normal0`, `normal1`, `normal2`: unit surface normal (out of surface, into domain)
 * - `momentum0`, `momentum1`, `momentum2` : momentum per volume
 * - `mass`: mass per volume (aka density)
 * - `energy`: total energy per volume
 * - `visc_stress0`, `visc_stress1`, `visc_stress2` : viscous stress at surface
 * - `mass_flux`: diffusive mass flux through surface
 * - `heat_flux`: surface heat flux
 */
void surface(Namespace&, Boundary_connection&, int i_fqpoint);

}
#endif
