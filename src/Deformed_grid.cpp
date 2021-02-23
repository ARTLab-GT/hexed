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
  n_vertices = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) n_vertices *= 2;
}

void Deformed_grid::add_vertices(std::vector<int> position, int i_dim)
{
  if (i_dim == n_dim)
  {
    vertices.emplace_back(vertices.size());
    Vertex& vertex = vertices.back();
    vertex.id_refs.push_back(vertex_ids.size());
    vertex_ids.push_back(vertex.id);
    vertex.parent_grid = this;
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
  for (int stride = 1; stride < n_vertices; stride *= 2)
  {
    for (int i_vertex = 0; i_vertex < n_vertices; ++i_vertex)
    {
      int i_neighbor = ((i_vertex/stride)%2 == 0) ? i_vertex + stride : i_vertex - stride;
      int id = vertex_ids[n_vertices*i_elem + i_vertex];
      vertices[id].neighbor_ids.push_back(n_vertices*i_elem + i_neighbor);
    }
  }
  for (int i = 0; i < n_qpoint/basis.rank*2*n_dim; ++i) node_adjustments.push_back(0);
  return i_elem;
}

}
