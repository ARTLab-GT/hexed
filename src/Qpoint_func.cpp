#include <Qpoint_func.hpp>
#include <Grid.hpp>
#include <math.hpp>

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

std::vector<double> Viscosity_func::operator()(Grid& grid, int i_element, int i_qpoint)
{
  Eigen::Map<Eigen::VectorXd> vert_visc (grid.element(i_element).viscosity(), custom_math::pow(2, grid.n_dim));
  Eigen::MatrixXd interp (grid.basis.row_size, 2);
  for (int i_node = 0; i_node < grid.basis.row_size; ++i_node) {
    interp(i_node, 0) = 1. - grid.basis.node(i_node);
    interp(i_node, 1) =      grid.basis.node(i_node);
  }
  Eigen::VectorXd qpoint_visc = custom_math::hypercube_matvec(interp, vert_visc);
  return {qpoint_visc[i_qpoint]};
}

}
