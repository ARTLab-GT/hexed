#ifndef CARTDG_QPOINT_POSITION_HPP_
#define CARTDG_QPOINT_POSITION_HPP_

#include <vector>
#include <Deformed_element.hpp>
#include <Basis.hpp>

namespace cartdg
{

/*
 * Returns a `std::vector` of size `n_dim` which contains the position of the `i_qpoint`th
 * quadrature point, based on the positions of the vertices. Warning: for Cartesian elements,
 * Jacobian will be the identity regardless of any deformity reflected in the qpoint position!
 */
std::array<double, 3> qpoint_position(Element&, const Basis&, int i_qpoint);

/*
 * Sets the Jacobian of a deformed element based on the current position of the quadrature points.
 * Jacobian is written to the `jacobian()` data member.
 * Also sets the `vertex_time_step_scale` member based on the Jacobian.
 */
void set_jacobian(Deformed_element&, const Basis&);

}
#endif
