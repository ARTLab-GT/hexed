#include <Qpoint_func.hpp>

namespace cartdg
{

std::vector<double> Jacobian_det_func::operator()(Element& element, int i_qpoint)
{
  return {element.jacobian_determinant(i_qpoint)};
}

}
