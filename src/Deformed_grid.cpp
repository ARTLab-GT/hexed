#include <Deformed_grid.hpp>

namespace cartdg
{

Deformed_grid::Deformed_grid(int n_var_arg, int n_dim_arg, int n_elem_arg,
                             double mesh_size_arg, Basis& basis_arg)
: Grid(n_var_arg, n_dim_arg, n_elem_arg, mesh_size_arg, basis_arg)
{
  if (n_elem_arg != 0)
  {
    auto message = "Capability to construct Deformed_grid with multiple elements is not implemented";
    throw message;
  }
}

void Deformed_grid::add_vertices(std::vector<int> position, int i_dim)
{
  if (i_dim == n_dim)
  {
    vertices.emplace_back(vertices.size());
    Vertex& vertex = vertices.back();
    vertex_ids.push_back(vertex.id);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      vertex.pos[i_dim] = position[i_dim]*mesh_size + origin[i_dim];
    }
  }
  else
  {
    add_vertices(position, i_dim + 1);
    ++position[i_dim];
    add_vertices(position, i_dim + 1);
  }
}

int Deformed_grid::add_element(std::vector<int> position)
{
  int i_elem = Grid::add_element(position);
  add_vertices(position, 0);
  for (int i = 0; i < n_qpoint/basis.rank*2*n_dim; ++i) node_adjustments.push_back(0);
  return i_elem;
}

}
