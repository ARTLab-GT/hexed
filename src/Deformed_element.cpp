#include <Deformed_element.hpp>
#include <math.hpp>

namespace cartdg
{

Deformed_element::Deformed_element(Storage_params params, std::vector<int> pos, double mesh_size)
: Element{params}, n_qpoint{params.n_qpoint()}, jac{n_dim*n_dim*n_qpoint}
{
  std::array<double, 3> first_pos;
  for (int i_dim = pos.size(); i_dim < 3; ++i_dim) pos.push_back(0);
  for (int i_dim = 0; i_dim < 3; ++i_dim) first_pos[i_dim] = pos[i_dim]*mesh_size;
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert)
  {
    auto vertex_pos = first_pos;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      int stride = custom_math::pow(2, n_dim - i_dim - 1);
      int i_row = (i_vert/stride)%2;
      vertex_pos[i_dim] += i_row*mesh_size;
    }
    vertices.emplace_back(vertex_pos);
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

}
