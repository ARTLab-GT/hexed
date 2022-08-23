#include <Vis_data.hpp>
#include <math.hpp>

namespace cartdg
{

Vis_data::Vis_data(Element& elem, const Qpoint_func&, const Basis& basis) :
  n_dim{elem.storage_params().n_dim},
  n_edge{custom_math::pow(2, n_dim - 1)*n_dim},
  row_size{elem.storage_params().row_size},
  n_qpoint{elem.storage_params().n_qpoint()},
  bas{basis},
  pos(n_qpoint*n_dim)
{
  // fetch data at quadrature points
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    auto p = elem.position(bas, i_qpoint);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) pos(n_qpoint*i_dim + i_qpoint) = p[i_dim];
  }
}

Eigen::MatrixXd Vis_data::edges(int n_div)
{
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_div + 1, 0., 1.))};
  const int nfqpoint = n_qpoint/row_size;
  Eigen::MatrixXd result(n_dim*(n_div + 1), n_edge);
  // interpolate qpoint position to edges
  Eigen::MatrixXd boundary {bas.boundary()};
  // interpolate all edges which point in direction `i_dim`
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    const int stride {custom_math::pow(row_size, n_dim - 1 - i_dim)};
    const int n_outer {n_qpoint/stride/row_size};
    // extract the `j_dim`th component of position
    for (int j_dim = 0; j_dim < n_dim; ++j_dim)
    {
      Eigen::MatrixXd edge_qpoints {row_size, n_edge/n_dim}; // quadrature points interpolated to the edge
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
        Eigen::VectorXd qpoint_slab {nfqpoint};
        for (int i_outer = 0; i_outer < n_outer; ++i_outer) {
          for (int i_inner = 0; i_inner < stride; ++i_inner) {
            qpoint_slab[i_outer*stride + i_inner] = pos(j_dim*n_qpoint + i_qpoint*stride + i_outer*stride*row_size + i_inner);
          }
        }
        // interpolate edge quadrature points to edge visualization points
        edge_qpoints.row(i_qpoint) = custom_math::hypercube_matvec(boundary, qpoint_slab);
      }
      result(Eigen::seqN(j_dim*(n_div + 1), n_div + 1),
             Eigen::seqN(i_dim*n_edge/n_dim, n_edge/n_dim)) = interp*edge_qpoints; // note: rhs auto resized to fit lhs
    }
  }
  return result;
}

}
