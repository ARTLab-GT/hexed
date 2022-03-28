#include <qpoint_position.hpp>
#include <math.hpp>

namespace cartdg
{

std::array<double, 3> qpoint_position(Element& element, Basis& basis, int i_qpoint)
{
  std::array<double, 3> pos;
  Storage_params par {element.storage_params()};
  const int n_vert = par.n_vertices();
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    Eigen::VectorXd vert_pos {n_vert};
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) vert_pos[i_vert] = element.vertex(i_vert).pos[i_dim];
    for (int j_dim = par.n_dim - 1; j_dim >= 0; --j_dim) {
      const int stride = custom_math::pow(par.row_size, par.n_dim - j_dim - 1);
      double node = basis.node((i_qpoint/stride)%par.row_size);
      Eigen::Matrix<double, 1, 2> interp {1. - node, node};
      vert_pos = custom_math::dimension_matvec(interp, vert_pos, j_dim);
    }
    pos[i_dim] = vert_pos[0];
  }
  return pos;
}

}
