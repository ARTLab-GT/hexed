#ifndef HEXED_HIL_PROPERTIES_HPP_
#define HEXED_HIL_PROPERTIES_HPP_

#include "Namespace.hpp"
#include "Element.hpp"

//! \brief functions that assign properties of objects to HIL variables
namespace hexed::hil_properties
{

/*! Assigns the follwing variables:
 * - `n_dim`: number of dimensions
 * - `is_def = elem.get_is_deformed()`
 * - `ref_level = elem.refinement_level()`
 * - `nom_sz` = elem.nominal_size()`
 * - `res_bad` = elem.resolution_badness()`
 * - `center0`, `center1`, `center2`: center of mass of vertices (not necessarily of the element itself)
 */
void element(Namespace&, Element& elem);

/*! Assigns the following variables `pos0`, `pos1`, `pos2`
 * to the position of the `i_qpoint`th quadrature point.
 * Trailing dimensions are set to 0.
 */
void position(Namespace&, Element&, const Basis&, int i_qpoint);

/*! Assigns the follwing variables:
 * - `state0`, `state1`, ... : conserved state variables
 * - `tss`: time step scale
 * - `art_visc`: artificial viscosity coefficient
 */
void state(Namespace&, Element&, int i_qpoint);

}
#endif
