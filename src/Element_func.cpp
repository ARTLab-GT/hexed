#include <Element_func.hpp>
#include <math.hpp>

namespace hexed
{

Elem_average::Elem_average(Qpoint_func& func) : qf{func} {}
Elem_l2::Elem_l2(Qpoint_func& func) : qf{func} {}

int Elem_average::n_var(int n_dim) const
{
  return qf.n_var(n_dim);
}

std::vector<double> avg(const Qpoint_func& qf, Element& elem, const Basis& basis, double time, int pow)
{
  auto params = elem.storage_params();
  auto weights = custom_math::pow_outer(basis.node_weights(), params.n_dim);
  int nv = qf.n_var(params.n_dim);
  std::vector<double> result(nv, 0.);
  double volume = 0.;
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
    auto qpoint = qf(elem, basis, i_qpoint, time);
    double w = weights[i_qpoint]*elem.jacobian_determinant(i_qpoint);
    for (int i_var = 0; i_var < nv; ++i_var) {
      result[i_var] += custom_math::pow(qpoint[i_var], pow)*w;
    }
    volume += w;
  }
  for (int i_var = 0; i_var < nv; ++i_var) result[i_var] /= volume;
  return result;
}

std::vector<double> Elem_average::operator()(Element& elem, const Basis& basis, double time) const
{
  return avg(qf, elem, basis, time, 1);
}

std::vector<double> Elem_l2::operator()(Element& elem, const Basis& basis, double time) const
{
  auto result = avg(qf, elem, basis, time, 2);
  for (double& val : result) val = std::sqrt(val);
  return result;
}

};
