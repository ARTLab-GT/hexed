#include <Qpoint_func.hpp>

namespace cartdg
{

Qpoint_func::Qpoint_func(Grid& g) : grid{g} {}

Jacobian_det_func::Jacobian_det_func(Grid& g) : Qpoint_func{g} {}

std::vector<double> Jacobian_det_func::operator()(int i_element, int i_qpoint)
{
  return {grid.element(i_element).jacobian_determinant(i_qpoint)};
}

}
