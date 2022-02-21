#include <Qpoint_func.hpp>
#include <Grid.hpp>

namespace cartdg
{

std::vector<double> Jacobian_det_func::operator()(Grid& grid, int i_element, int i_qpoint)
{
  return {grid.element(i_element).jacobian_determinant(i_qpoint)};
}

std::vector<double> Time_step_scale_func::operator()(Grid& grid, int i_element, int i_qpoint)
{
  return {grid.element(i_element).time_step_scale()[i_qpoint]};
}

}
