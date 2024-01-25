#include <Qpoint_func.hpp>
#include <math.hpp>
#include <utils.hpp>
#include <Element.hpp>
#include <hil_properties.hpp>

namespace hexed
{

Qpoint_expr::Qpoint_expr(Struct_expr expr, const Interpreter& inter) : _expr{expr}, _inter{inter} {}

std::vector<double> Qpoint_expr::operator()(Element& elem, const Basis& basis, int i_qpoint, double time) const
{
  auto sub = _inter.make_sub();
  hil_properties::element(*sub.variables, elem);
  hil_properties::position(*sub.variables, elem, basis, i_qpoint);
  hil_properties::state(*sub.variables, elem, i_qpoint);
  sub.variables->assign("time", time);
  return _expr.eval(sub);
}

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

std::vector<double> Physical_residual::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  std::vector<double> result;
  const int n_qpoint = element.storage_params().n_qpoint();
  for (int i_var = 0; i_var < element.storage_params().n_var; ++i_var) {
    result.push_back(element.residual_cache()[i_var*n_qpoint + i_qpoint]/element.time_step_scale()[i_qpoint]);
  }
  return result;
}

std::vector<double> Art_visc_coef::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.art_visc_coef()[i_qpoint]};
}

std::vector<double> Fix_admis_coef::operator()(Element& element, const Basis&, int i_qpoint, double time) const
{
  return {element.fix_admis_coef()[i_qpoint]};
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
  for (double& r : result) r = math::pow(r, exp);
  return result;
}

std::vector<double> Advection_state::operator()(Element& elem, const Basis&, int i_qpoint, double time) const
{
  int elem_rs = elem.storage_params().row_size;
  if (elem_rs != rs) {
    throw std::runtime_error(format_str(100, "`Advection_state` initialized with row size %d but called with row size %d", rs, elem_rs));
  }
  std::vector<double> as;
  for (int i_as = 0; i_as < rs; ++i_as) {
    as.push_back(elem.advection_state()[i_as*elem.storage_params().n_qpoint() + i_qpoint]);
  }
  return as;
}

int Art_visc_forcing::n_var(int n_dim) const {return 4;}

std::vector<double> Art_visc_forcing::operator()(Element& elem, const Basis&, int i_qpoint, double time) const
{
  std::vector<double> forcing;
  for (int i_forcing = 0; i_forcing < elem.storage_params().n_forcing; ++i_forcing) {
    forcing.push_back(elem.art_visc_forcing()[i_forcing*elem.storage_params().n_qpoint() + i_qpoint]);
  }
  return forcing;
}

}
