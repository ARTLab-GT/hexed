#include <Qpoint_func.hpp>
#include <Grid.hpp>
#include <math.hpp>

namespace cartdg
{

std::vector<double> Qpoint_func::operator()(Grid& grid, int i_element, int i_qpoint) const
{
  return operator()(grid.element(i_element), grid.basis, i_qpoint, grid.time);
}

std::vector<double> Jacobian_det_func::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.jacobian_determinant(i_qpoint)};
}

std::vector<double> Time_step_scale_func::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.time_step_scale()[i_qpoint]};
}

std::vector<double> Viscosity_func::operator()(Element& element, const Basis& basis, int i_qpoint, double time) const
{
  Storage_params params {element.storage_params()};
  Eigen::Map<Eigen::VectorXd> vert_visc (element.viscosity(), custom_math::pow(2, params.n_dim));
  Eigen::MatrixXd interp (params.row_size, 2);
  for (int i_node = 0; i_node < basis.row_size; ++i_node) {
    interp(i_node, 0) = 1. - basis.node(i_node);
    interp(i_node, 1) =      basis.node(i_node);
  }
  Eigen::VectorXd qpoint_visc = custom_math::hypercube_matvec(interp, vert_visc);
  return {qpoint_visc[i_qpoint]};
}

std::vector<double> Physical_update::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  std::vector<double> result;
  const int n_qpoint = element.storage_params().n_qpoint();
  for (int i_var = 0; i_var < element.storage_params().n_var; ++i_var) {
    result.push_back((element.stage(0)[i_var*n_qpoint + i_qpoint] - element.stage(1)[i_var*n_qpoint + i_qpoint])/element.time_step_scale()[i_qpoint]);
  }
  return result;
}

}
