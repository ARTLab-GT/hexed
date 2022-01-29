#include <Deformed_element.hpp>
#include <math.hpp>

namespace cartdg
{

Deformed_element::Deformed_element(Storage_params params, std::vector<int> pos, double mesh_size)
: Element{params}, n_qpoint{params.n_qpoint()}, jac{n_dim*n_dim*n_qpoint},
  node_adj{Eigen::VectorXd::Zero(n_qpoint/params.row_size*n_dim*2)}
{
  // set position of vertex 0
  std::array<double, 3> first_pos;
  int n_pos_set = std::min<int>(pos.size(), n_dim);
  for (int i_dim = 0; i_dim < n_pos_set; ++i_dim) first_pos[i_dim] = pos[i_dim]*mesh_size;
  for (int i_dim = n_pos_set; i_dim < 3; ++i_dim) first_pos[i_dim] = 0.;
  // construct vertices
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert)
  {
    // compute position of vertex
    auto vertex_pos = first_pos;
    int stride [3];
    int i_row [3];
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      stride[i_dim] = custom_math::pow(2, n_dim - i_dim - 1);
      i_row[i_dim] = (i_vert/stride[i_dim])%2;
      vertex_pos[i_dim] += i_row[i_dim]*mesh_size;
    }
    vertices.emplace_back(vertex_pos);
    // establish vertex connections (that is, edges).
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      if (i_row[i_dim]) Vertex::connect(*vertices.back(), *vertices[i_vert - stride[i_dim]]);
    }
  }
}

double* Deformed_element::jacobian()
{
  return jac.data();
}

double Deformed_element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return jac[(n_dim*i_dim + j_dim)*n_qpoint + i_qpoint];
}

double* Deformed_element::node_adjustments()
{
  return node_adj.data();
}

Deformed_face::Deformed_face(Storage_params params)
: n_dim(params.n_dim), n_fqpoint(params.n_qpoint()/params.row_size), jac(n_dim*n_dim*n_fqpoint)
{
}

double* Deformed_face::jacobian()
{
  return jac.data();
}

double Deformed_face::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return jac((i_dim*n_dim + j_dim)*n_fqpoint + i_qpoint);
}

}
