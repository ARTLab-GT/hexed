#ifndef HEXED_VIS_VARIABLES_HPP_
#define HEXED_VIS_VARIABLES_HPP_

#include "Namespace.hpp"
#include "Element.hpp"
#include "connection.hpp"

//! \brief functions that assign visualization data to HIL variables
namespace hexed::vis_variables {

/*! \brief Assigns the variables `pos0`, `pos1`, `pos2`
 * \details to the position of the `i_qpoint`th quadrature point.
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
 * \attention If `Solver::update` or `Solver::update_art_visc_smoothness` have been called since the last call to `Solver::compute_residuals`
 * then the residual variables will be incorrect.
 */
void state(Namespace&, Element&);

}
#endif
