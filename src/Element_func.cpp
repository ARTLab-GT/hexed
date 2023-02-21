#include <Element_func.hpp>
#include <math.hpp>

namespace hexed
{

std::vector<double> Element_func::operator()(Element& elem, const Basis& basis, int, double time) const
{
  return operator()(elem, basis, time);
}

Elem_average::Elem_average    (const Qpoint_func& func) : qf{func} {}
Elem_l2::Elem_l2              (const Qpoint_func& func) : qf{func} {}
Elem_nonsmooth::Elem_nonsmooth(const Qpoint_func& func) : qf{func} {}

int Elem_average::n_var(int n_dim) const
{
  return qf.n_var(n_dim);
}

//! \private computes the average of some function raised to a power
std::vector<double> avg(const Qpoint_func& qf, Element& elem, const Basis& basis, double time, int pow)
{
  auto params = elem.storage_params();
  auto weights = math::pow_outer(basis.node_weights(), params.n_dim);
  int nv = qf.n_var(params.n_dim);
  std::vector<double> result(nv, 0.);
  double volume = 0.;
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
    auto qpoint = qf(elem, basis, i_qpoint, time);
    double w = weights[i_qpoint]*elem.jacobian_determinant(i_qpoint);
    for (int i_var = 0; i_var < nv; ++i_var) {
      result[i_var] += math::pow(qpoint[i_var], pow)*w;
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

std::vector<double> Elem_nonsmooth::operator()(Element& elem, const Basis& basis, double time) const
{
  auto params = elem.storage_params();
  auto weights = math::pow_outer(basis.node_weights(), params.n_dim - 1);
  const int nv = qf.n_var(params.n_dim);
  const int nq = params.n_qpoint();
  std::vector<double> result(nv, 0.);
  // find the value of the function we're taking the nonsmoothness of
  Eigen::MatrixXd vals(nq, nv);
  for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
    auto qpoint = qf(elem, basis, i_qpoint, time);
    for (int i_var = 0; i_var < nv; ++i_var) vals(i_qpoint, i_var) = qpoint[i_var];
  }
  // projector onto the highest-order hierarchical orthogonal basis function
  Eigen::MatrixXd max_orth = (basis.orthogonal(params.row_size - 1).cwiseProduct(basis.node_weights())).transpose();
  for (int i_var = 0; i_var < nv; ++i_var) {
    // compute RMS of the following for each dimension:
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
      // take the highest-order component in each dimension...
      Eigen::VectorXd proj = math::dimension_matvec(max_orth, vals(Eigen::all, i_var), i_dim);
      // ...and compute L2 norm of it in the other dimensions
      result[i_var] += proj.dot(proj.cwiseProduct(weights));
    }
    result[i_var] = std::sqrt(result[i_var]/params.n_dim);
  }
  return result;
}

};
