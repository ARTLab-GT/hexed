#include <Deformed_grid.hpp>
#include <get_local_deformed_cpg_euler.hpp>
#include <get_neighbor_deformed_cpg_euler.hpp>
#include <get_gbc_cpg_euler.hpp>
#include <get_nonpen_cpg_euler.hpp>

namespace cartdg
{

Deformed_grid::Deformed_grid(int n_var_arg, int n_dim_arg, int n_elem_arg,
                             double mesh_size_arg, Basis& basis_arg)
: Grid(n_var_arg, n_dim_arg, n_elem_arg, mesh_size_arg, basis_arg)
{
  if (n_elem_arg != 0)
  {
    throw std::runtime_error("Capability to construct Deformed_grid with multiple elements is not implemented.");
  }
  n_vertices = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) n_vertices *= 2;
  neighbor_storage.resize(3);

  default_jacobian.clear();
  default_jacobian.resize(n_dim*n_dim*n_qpoint, 0.);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      default_jacobian[i_dim*(n_dim + 1)*n_qpoint + i_qpoint] = 1.;
    }
  }
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

void Deformed_grid::add_wall(int i_elem, int i_dim, bool is_positive_face)
{
  i_elem_wall.push_back(i_elem);
  i_dim_wall.push_back(i_dim);
  is_positive_wall.push_back((int)is_positive_face);
}

double Deformed_grid::jacobian_det(int i_elem, int i_qpoint)
{
  Eigen::MatrixXd jac_mat (n_dim, n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim)
    {
      int i_value = ((i_elem*n_dim + i_dim)*n_dim + j_dim)*n_qpoint + i_qpoint;
      jac_mat(i_dim, j_dim) = jacobian[i_value];
    }
  }
  return jac_mat.fullPivLu().determinant();
}

void Deformed_grid::execute_local(Kernel_settings& settings)
{
  get_local_deformed_cpg_euler(n_dim, basis.rank)(state_r(), state_w(), jacobian.data(), n_elem,
                                                  basis, settings);
}

void Deformed_grid::execute_neighbor(Kernel_settings& settings)
{
  get_neighbor_deformed_cpg_euler(n_dim, basis.rank)(state_connections_r(), state_connections_w(), jacobian_neighbors.data(), neighbor_axes.data(), neighbor_is_positive.data(), neighbor_storage[0].size()/2, basis, settings);
  get_gbc_cpg_euler(n_dim, basis.rank)(ghost_bound_conds, state_r(), state_w(), basis, settings);
  get_nonpen_cpg_euler(n_dim, basis.rank)(state_r(), state_w(), jacobian.data(), i_elem_wall.data(), i_dim_wall.data(), is_positive_wall.data(), i_elem_wall.size(), basis, settings);
}

void Deformed_grid::execute_req_visc(Kernel_settings& settings)
{
}

void Deformed_grid::execute_cont_visc(Kernel_settings& settings)
{
}

void Deformed_grid::execute_local_derivative(int i_var, int i_axis, Kernel_settings& settings)
{
  double* sr = state_r();
  double* sw = state_w();
  for (int i_data = 0; i_data < n_dof*n_elem; ++i_data)
  {
    sw[i_data] = sr[i_data];
  }
}

void Deformed_grid::execute_neighbor_derivative(int i_var, int i_axis, Kernel_settings& settings)
{
}

void Deformed_grid::execute_av_flux(Kernel_settings& settings)
{
}

void Deformed_grid::execute_local_av(int i_var, int i_axis, Kernel_settings& settings)
{
}

void Deformed_grid::execute_neighbor_av(int i_var, int i_axis, Kernel_settings& settings)
{
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
  for (int i_side : {0, 1})
  {
    neighbor_inds.push_back(i_elem[i_side]);
    neighbor_axes.push_back(i_axis[i_side]);
    neighbor_is_positive.push_back(is_positive[i_side]);
  }
}

void Deformed_grid::update_connections()
{
  for (int i_elem : neighbor_inds)
  {
    jacobian_neighbors.push_back(jacobian.data() + i_elem*n_dim*n_dim*n_qpoint);
    for (int i_stage = 0; i_stage < 3; ++i_stage)
    {
      double* state_data = state_storage[i_stage].data();
      neighbor_storage[i_stage].push_back(state_data + n_dof*i_elem);
    }
  }
}

void Deformed_grid::connect_non_def(std::array<int, 2> i_elem, std::array<int, 2> i_axis,
                                    std::array<bool, 2> is_positive, Grid& other_grid)
{
  for (int i_stage = 0; i_stage < 3; ++i_stage)
  {
    neighbor_storage[i_stage].push_back(state_storage[i_stage].data() + n_dof*i_elem[0]);
    neighbor_storage[i_stage].push_back(other_grid.state_storage[i_stage].data()
                                        + other_grid.n_dof*i_elem[1]);
  }
  for (int i_side : {0, 1})
  {
    neighbor_axes.push_back(i_axis[i_side]);
    neighbor_is_positive.push_back(int(is_positive[i_side]));
  }
  jacobian_neighbors.push_back(jacobian.data() + i_elem[0]*n_dim*n_dim*n_qpoint);
  jacobian_neighbors.push_back(default_jacobian.data());
}

void Deformed_grid::visualize(std::string file_name)
{
  Grid::visualize(file_name + "_deformed");
}

std::vector<double> Deformed_grid::face_integral(Domain_func& integrand, int i_elem, int i_dim, bool is_positive)
{
  auto pos = get_pos(i_elem);
  double* elem_state = state_r() + i_elem*n_dof;
  auto row_weights = basis.node_weights();
  int stride = std::pow(basis.rank, n_dim - 1 - i_dim);
  std::vector<double> total;
  for (int i_outer = 0; i_outer < n_qpoint/basis.rank/stride; ++i_outer)
  {
    for (int i_inner = 0; i_inner < stride; ++i_inner)
    {
      int i_qpoint = (i_outer*basis.rank + int(is_positive)*(basis.rank - 1.))*stride + i_inner;
      std::vector<double> qpoint_pos;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        qpoint_pos.push_back(pos[i_qpoint + n_qpoint*j_dim]);
      }
      std::vector<double> qpoint_state;
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        qpoint_state.push_back(elem_state[i_qpoint + i_var*n_qpoint]);
      }
      auto qpoint_integrand = integrand(qpoint_pos, time, qpoint_state);
      if (total.size() < qpoint_integrand.size())
      {
        total.resize(qpoint_integrand.size());
      }
      Eigen::MatrixXd qpoint_jacobian (n_dim, n_dim);
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        for (int k_dim = 0; k_dim < n_dim; ++k_dim)
        {
          qpoint_jacobian(j_dim, k_dim) = jacobian[((i_elem*n_dim + j_dim)*n_dim + k_dim)*n_qpoint + i_qpoint];
        }
        qpoint_jacobian(j_dim, i_dim) = 0.;
      }
      double face_jac_det = 0.;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        qpoint_jacobian(j_dim, i_dim) = 1.;
        double dim_jac = qpoint_jacobian.determinant();
        face_jac_det += dim_jac*dim_jac;
        qpoint_jacobian(j_dim, i_dim) = 0.;
      }
      face_jac_det = std::sqrt(face_jac_det);
      int i_face = i_inner + i_outer*stride;
      double weight = 1.;
      for (int j_dim = 0; j_dim < n_dim - 1; ++j_dim)
      {
        int stride_j = std::pow(basis.rank, j_dim);
        weight *= row_weights((i_face/stride_j)%basis.rank);
      }
      for (int i_var = 0; i_var < int(total.size()); ++i_var)
      {
        total[i_var] += weight*face_jac_det*qpoint_integrand[i_var];
      }
    }
  }
  for (int i_var = 0; i_var < int(total.size()); ++i_var)
  {
    total[i_var] *= std::pow(mesh_size, n_dim - 1);
  }
  return total;
}

std::vector<double> Deformed_grid::surface_integral(Domain_func& integrand)
{
  std::vector<double> total;
  for (int i_wall = 0; i_wall < int(i_elem_wall.size()); ++i_wall)
  {
    auto fi = face_integral(integrand, i_elem_wall[i_wall], i_dim_wall[i_wall], is_positive_wall[i_wall]);
    if (total.size() < fi.size())
    {
      total.resize(fi.size());
    }
    for (int i_var = 0; i_var < int(total.size()); ++i_var)
    {
      total[i_var] += fi[i_var];
    }
  }
  return total;
}

void Deformed_grid::calc_jacobian()
{
  jacobian.resize(n_elem*n_dim*n_dim*n_qpoint);
  // FIXME: make this a kernel
  auto diff_mat = basis.diff_mat();
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> elem_pos = get_pos(i_elem);
    for (int i_dim = 0, stride = n_qpoint/basis.rank; i_dim < n_dim; ++i_dim, stride /= basis.rank)
    {
      for (int i_outer = 0; i_outer < n_qpoint/(stride*basis.rank); ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          for (int j_dim = 0; j_dim < n_dim; ++j_dim)
          {
            Eigen::VectorXd row_pos (basis.rank);
            int row_start = i_outer*stride*basis.rank + i_inner;
            for (int i_qpoint = 0; i_qpoint < basis.rank; ++i_qpoint)
            {
              row_pos(i_qpoint) = elem_pos[j_dim*n_qpoint + row_start + i_qpoint*stride];
            }
            auto row_jacobian = diff_mat*row_pos;
            for (int i_qpoint = 0; i_qpoint < basis.rank; ++i_qpoint)
            {
              jacobian[((i_elem*n_dim + j_dim)*n_dim + i_dim)*n_qpoint
                       + row_start + i_qpoint*stride] = row_jacobian(i_qpoint)/mesh_size;
            }
          }
        }
      }
    }
  }
}

double** Deformed_grid::state_connections_r()
{
  return neighbor_storage[i_read].data();
}

double** Deformed_grid::state_connections_w()
{
  return neighbor_storage[i_write].data();
}

}
