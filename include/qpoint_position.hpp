#ifndef CARTDG_QPOINT_POSITION_HPP_
#define CARTDG_QPOINT_POSITION_HPP_

#include <vector>
#include <Element.hpp>
#include <Basis.hpp>

namespace cartdg
{

/*
 * Returns a `std::vector` of size `n_dim` which contains the position of the `i_qpoint`th
 * quadrature point, based on the positions of the vertices. Warning: for Cartesian elements,
 * Jacobian will be the identity regardless of any deformity reflected in the qpoint position!
 */
std::array<double, 3> qpoint_position(const Element&, const Basis&, int i_qpoint);

}
#endif
