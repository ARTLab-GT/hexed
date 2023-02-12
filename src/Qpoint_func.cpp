#include <Qpoint_func.hpp>
#include <math.hpp>

namespace hexed
{

std::vector<double> Jacobian_det_func::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.jacobian_determinant(i_qpoint)};
}

std::vector<double> Jac_inv_det_func::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {1./element.jacobian_determinant(i_qpoint)};
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

std::vector<double> Art_visc_coef::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.art_visc_coef()[i_qpoint]};
}

std::string Pow::variable_name(int n_dim, int i_var) const
{
  char buffer [1000];
  snprintf(buffer, 1000, "(%s)^%i", qf.variable_name(n_dim, i_var).c_str(), exp);
  return buffer;
}

std::vector<double> Pow::operator()(Element& e, const Basis& b, int i_qpoint, double time) const
{
  std::vector<double> result = qf(e, b, i_qpoint, time);
  for (double& r : result) r = custom_math::pow(r, exp);
  return result;
}

int Qf_concat::n_var(int n_dim) const
{
  return 0;
}

std::string Qf_concat::variable_name(int n_dim, int i_var) const
{
  return "";
}

std::vector<double> Qf_concat::operator()(Element&, const Basis&, int i_qpoint, double time) const
{
  return {};
}

}
