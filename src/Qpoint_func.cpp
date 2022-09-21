#include <Qpoint_func.hpp>
#include <math.hpp>

namespace hexed
{

std::vector<double> Jacobian_det_func::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.jacobian_determinant(i_qpoint)};
}

std::vector<double> Time_step_scale_func::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.time_step_scale()[i_qpoint]};
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
