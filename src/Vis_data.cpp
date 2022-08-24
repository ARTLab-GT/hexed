#include <Vis_data.hpp>
#include <math.hpp>

namespace cartdg
{

Vis_data::Vis_data(Element& elem, const Qpoint_func& func, const Basis& basis, double time) :
  n_dim{elem.storage_params().n_dim},
  n_edge{custom_math::pow(2, n_dim - 1)*n_dim},
  row_size{elem.storage_params().row_size},
  n_qpoint{elem.storage_params().n_qpoint()},
  n_var{func.n_var(n_dim)},
  bas{basis},
  vars(n_qpoint*n_var)
{
  // fetch data at quadrature points
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    auto v = func(elem, bas, i_qpoint, time);
    for (int i_var = 0; i_var < n_var; ++i_var) vars(n_qpoint*i_var + i_qpoint) = v[i_var];
  }
}

Eigen::VectorXd Vis_data::edges(int n_sample)
{
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  const int nfqpoint = n_qpoint/row_size;
  Eigen::MatrixXd result(n_var*n_sample, n_edge);
  // interpolate qpoint varsition to edges
  Eigen::MatrixXd boundary {bas.boundary()};
  // interpolate all edges which point in direction `i_dim`
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    const int stride {custom_math::pow(row_size, n_dim - 1 - i_dim)};
    const int n_outer {n_qpoint/stride/row_size};
    // extract the `i_var`th variable
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      Eigen::MatrixXd edge_qpoints {row_size, n_edge/n_dim}; // quadrature points interpolated to the edge
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
        Eigen::VectorXd qpoint_slab {nfqpoint};
        for (int i_outer = 0; i_outer < n_outer; ++i_outer) {
          for (int i_inner = 0; i_inner < stride; ++i_inner) {
            qpoint_slab[i_outer*stride + i_inner] = vars(i_var*n_qpoint + i_qpoint*stride + i_outer*stride*row_size + i_inner);
          }
        }
        // interpolate edge quadrature points to edge visualization points
        edge_qpoints.row(i_qpoint) = custom_math::hypercube_matvec(boundary, qpoint_slab);
      }
      result(Eigen::seqN(i_var*n_sample, n_sample),
             Eigen::seqN(i_dim*n_edge/n_dim, n_edge/n_dim)) = interp*edge_qpoints; // note: rhs auto resized to fit lhs
    }
  }
  result.resize(result.size(), 1);
  return result;
}

Eigen::VectorXd Vis_data::interior(int n_sample)
{
  const int n_block = custom_math::pow(n_sample, n_dim);
  return Eigen::VectorXd::Zero(n_block*n_var);
}

}
