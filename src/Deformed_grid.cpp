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

Vertex& Deformed_grid::get_vertex(int i_vertex)
{
  return vertices[vertex_ids[i_vertex]];
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

std::vector<double> Deformed_grid::get_pos(int i_elem)
{
  std::vector<double> elem_pos (n_qpoint*n_dim);
  for (int i_vertex = 0; i_vertex < n_vertices; ++i_vertex)
  {
    int id = vertex_ids[n_vertices*i_elem + i_vertex];
    int i_node = 0;
    for (int vertex_stride = 1, node_stride = 1;
         vertex_stride < n_vertices;
         vertex_stride *= 2, node_stride *= basis.rank)
    {
      if ((i_vertex/vertex_stride)%2 == 1) i_node += node_stride*(basis.rank - 1);
    }
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      elem_pos[n_qpoint*i_dim + i_node] = get_vertex(id).pos[i_dim];
    }
  }

  for (int i_dim = 0, stride = n_qpoint/basis.rank; i_dim < n_dim; ++i_dim, stride /= basis.rank)
  {
    for (int i_node = 0; i_node < n_qpoint; ++i_node)
    {
      int coord = (i_node/stride)%basis.rank;
      double dist = basis.node(coord);
      int i_node0 = i_node - coord*stride;
      int i_node1 = i_node0 + (basis.rank - 1)*stride;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        int i = j_dim*n_qpoint;
        elem_pos[i + i_node] = (1. - dist)*elem_pos[i + i_node0] + dist*elem_pos[i + i_node1];
      }
    }
  }

  std::vector<double> warped_elem_pos = elem_pos;
  for (int i_dim = 0, stride = n_qpoint/basis.rank; i_dim < n_dim; ++i_dim, stride /= basis.rank)
  {
    for (int i_node = 0; i_node < n_qpoint; ++i_node)
    {
      int coord = (i_node/stride)%basis.rank;
      int i_node0 = i_node - coord*stride;
      int i_node1 = i_node0 + (basis.rank - 1)*stride;
      int i_adjust = 2*(i_dim + i_elem*n_dim)*n_qpoint/basis.rank
                     + i_node/(stride*basis.rank)*stride + i_node%stride;
      double adjust0 = node_adjustments[i_adjust];
      double adjust1 = node_adjustments[i_adjust + n_qpoint/basis.rank];
      double dist = basis.node(coord);
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        int i = j_dim*n_qpoint;
        warped_elem_pos[i + i_node] += (elem_pos[i + i_node1] - elem_pos[i + i_node0])*((1. - dist)*adjust0 + dist*adjust1);
      }
    }
  }
  return warped_elem_pos;
}

void Deformed_grid::connect(std::array<int, 2> i_elem, std::array<int, 2> i_axis,
                            std::array<bool, 2> is_positive)
{
  std::array<std::vector<int>, 2> id_inds;
  std::array<int, 2> strides;
  for (int i_side : {0, 1})
  {
    int stride = n_vertices/2;
    for (int i = 0; i < i_axis[i_side]; ++i) stride /= 2;
    strides[i_side] = stride;
    for (int i_vertex = 0; i_vertex < n_vertices; ++i_vertex)
    {
      if ((i_vertex/stride)%2 == int(is_positive[i_side]))
      {
        id_inds[i_side].push_back(i_vertex + i_elem[i_side]*n_vertices);
      }
    }
  }
  if ((is_positive[0] != is_positive[1]) && (i_axis[0] != i_axis[1]))
  {
    if (n_dim == 3)
    {
      int stride = strides[0];
      if (i_axis[0] < i_axis[1]) stride /= 2;
      for (int i : {0, 1}) std::swap(id_inds[1][i*2/stride], id_inds[1][i*2/stride + stride]);
    }
    else std::swap(id_inds[1][0], id_inds[1][1]);
  }
  if ((i_axis[0] == 0 && i_axis[1] == 2) || (i_axis[0] == 2 && i_axis[1] == 0))
  {
    std::swap(id_inds[1][1], id_inds[1][2]);
  }
  for (int i_vertex = 0; i_vertex < n_vertices/2; ++i_vertex)
  {
    get_vertex(id_inds[0][i_vertex]).eat(get_vertex(id_inds[1][i_vertex]));
  }
}

void Deformed_grid::visualize(std::string file_name)
{
  Grid::visualize(file_name + "_deformed");
}

}
