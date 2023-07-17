#include <Element_func.hpp>
#include <math.hpp>
#include <connection.hpp>

namespace hexed
{

std::vector<double> Element_func::operator()(Element& elem, const Basis& basis, int, double time) const
{
  return operator()(elem, basis, time);
}

std::vector<double> Element_info::operator()(Element& elem, const Basis&, double time) const
{
  return operator()(elem);
}

std::vector<double> Element_info::operator()(Boundary_connection& con, int i_fqpoint, double time) const
{
  return operator()(con.element());
}

std::vector<double> Resolution_badness::operator()(Element& elem) const {return {elem.resolution_badness};}

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

std::vector<double> Equiangle_skewness::operator()(Element& elem, const Basis& basis, double time) const
{
  const int nq = elem.storage_params().n_qpoint();
  const int nd = elem.storage_params().n_dim;
  const int nv = elem.storage_params().n_vertices();
  // fetch Jacobian at quadrature points
  Mat<dyn, dyn> qpoint_jac(nq, nd*nd);
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    for (int j_dim = 0; j_dim < nd; ++j_dim) {
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        qpoint_jac(i_qpoint, i_dim*nd + j_dim) = elem.jacobian(i_dim, j_dim, i_qpoint);
      }
    }
  }
  // extrapolate Jacobian to vertices
  Mat<dyn, dyn> vertex_jac(nv, nd*nd);
  for (int i_jac = 0; i_jac < nd*nd; ++i_jac) {
    vertex_jac(all, i_jac) = math::hypercube_matvec(basis.boundary(), qpoint_jac(all, i_jac));
  }
  // evaluate skewness at vertices
  double max = 0.;
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    // fetch vertex jacobian
    Mat<dyn, dyn> jac(nd, nd);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      for (int j_dim = 0; j_dim < nd; ++j_dim) {
        jac(i_dim, j_dim) = vertex_jac(i_vert, i_dim*nd +j_dim);
      }
    }
    // for every combination of dimensions...
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      for (int j_dim = i_dim + 1; j_dim < nd; ++j_dim) {
        // evaluate equiangle criterion
        double angle = std::acos(jac(all, i_dim).normalized().dot(jac(all, j_dim).normalized()));
        max = std::max(max, std::abs(2*angle/M_PI - 1));
      }
    }
  }
  return {max};
}

std::vector<double> Is_deformed::operator()(Element& elem) const {return {double(elem.get_is_deformed())};}
std::vector<double> Has_tree::operator()(Element& elem) const {return {double(bool(elem.tree))};}

};
