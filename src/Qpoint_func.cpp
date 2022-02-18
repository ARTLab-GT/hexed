#include <Qpoint_func.hpp>

namespace cartdg
{

std::vector<double> Jacobian_det_func::operator()(Grid& grid, int i_element, int i_qpoint)
{
  return {grid.element(i_element).jacobian_determinant(i_qpoint)};
}

}
