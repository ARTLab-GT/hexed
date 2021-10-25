#include <Deformed_grid.hpp>
#include <Regular_grid.hpp>
#include <get_mcs_deformed_convective.hpp>
#include <get_write_face_deformed.hpp>
#include <get_neighbor_deformed_convective.hpp>
#include <get_gbc_convective.hpp>
#include <get_local_deformed_convective.hpp>
#include <get_nonpen_convective.hpp>
#include <get_req_visc_deformed_convective.hpp>
#include <math.hpp>

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
}

Element& Deformed_grid::element(int i_elem)
{
  return *elements[i_elem];
}

Deformed_element& Deformed_grid::deformed_element(int i_elem)
{
  return *elements[i_elem];
}

double Deformed_grid::stable_time_step(double cfl_by_stable_cfl, Kernel_settings& settings)
{
  double cfl = cfl_by_stable_cfl*get_stable_cfl();
  settings.i_read = i_read;
  return cfl*mesh_size/get_mcs_deformed_convective(n_dim, basis.row_size)(elements, settings);
}

Deformed_elem_con Deformed_grid::connection(int i_con)
{
  return elem_cons[i_con];
}

Deformed_elem_wall Deformed_grid::def_elem_wall(int i_wall)
{
  return walls[i_wall];
}

int Deformed_grid::add_element(std::vector<int> position)
{
  int i_elem = Grid::add_element(position);
  elements.emplace_back(new Deformed_element {storage_params, position, mesh_size});
  for (int i_vert = 0; i_vert < n_vertices; ++i_vert)
  {
    Vertex& vert = elements.back()->vertex(i_vert);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) vert.pos[i_dim] += origin[i_dim];
    Vertex::Non_transferable_ptr ptr {vert};
    vertices.push_back(vert);
  }
  return i_elem;
}

std::vector<double> Deformed_grid::get_pos(int i_elem)
{
  std::vector<double> elem_pos (n_qpoint*n_dim);
  Deformed_element& elem {deformed_element(i_elem)};
  Eigen::MatrixXd lin_interp {basis.row_size, 2};
  for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint)
  {
    double node {basis.node(i_qpoint)};
    lin_interp(i_qpoint, 0) = 1. - node;
    lin_interp(i_qpoint, 1) = node;
  }

  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    Eigen::VectorXd vert_pos {n_vertices};
    for (int i_vertex = 0; i_vertex < n_vertices; ++i_vertex)
    {
      vert_pos[i_vertex] = elem.vertex(i_vertex).pos[i_dim];
    }
    Eigen::Map<Eigen::VectorXd> dim_pos {elem_pos.data() + i_dim*n_qpoint, n_qpoint};
    dim_pos = custom_math::hypercube_matvec(lin_interp, vert_pos);
  }

  std::vector<double> warped_elem_pos = elem_pos;
  Eigen::MatrixXd boundary {basis.boundary()};
  for (int i_dim = 0, stride = n_qpoint/basis.row_size; i_dim < n_dim; ++i_dim, stride /= basis.row_size)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint/basis.row_size; ++i_qpoint)
    {
      int i_row_start = i_qpoint/stride*stride*basis.row_size + i_qpoint%stride;
      double* node_adj = elem.node_adjustments();
      int i_adjust = 2*i_dim*n_qpoint/basis.row_size + i_qpoint;
      Eigen::VectorXd adjust {{node_adj[i_adjust], node_adj[i_adjust + n_qpoint/basis.row_size]}};
      typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> stride_t;
      typedef Eigen::Map<Eigen::MatrixXd, 0, stride_t> map;
      map orig {&elem_pos[i_row_start], basis.row_size, n_dim, stride_t{n_qpoint, stride}};
      map warped {&warped_elem_pos[i_row_start], basis.row_size, n_dim, stride_t{n_qpoint, stride}};
      Eigen::MatrixXd face_pos = boundary*orig;
      warped += lin_interp*(adjust*(face_pos.row(1) - face_pos.row(0)));
    }
  }
  return warped_elem_pos;
}

void Deformed_grid::add_wall(int i_elem, int i_dim, bool is_positive_face)
{
  Deformed_elem_wall wall {elements[i_elem].get(), i_dim, is_positive_face, i_elem};
  walls.push_back(wall);
}

void Deformed_grid::execute_write_face(Kernel_settings& settings)
{
  settings.i_read = i_read;
  get_write_face_deformed(n_dim, basis.row_size)(elements, basis, settings);
}

void Deformed_grid::execute_neighbor(Kernel_settings& settings)
{
  get_neighbor_deformed_convective(n_dim, basis.row_size)(elem_cons, basis, settings);
  get_gbc_convective(n_dim, basis.row_size)(*this, basis, settings);
  get_nonpen_convective(n_dim, basis.row_size)(walls, basis, settings);
}

void Deformed_grid::execute_local(Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_local_deformed_convective(n_dim, basis.row_size)(elements, basis, settings);
}

void Deformed_grid::execute_req_visc(Kernel_settings& settings)
{
  get_req_visc_deformed_convective(n_dim, basis.row_size)(elements, basis, settings);
}

void Deformed_grid::execute_cont_visc(Kernel_settings& settings)
{
}

void Deformed_grid::execute_local_derivative(int i_var, int i_dim, Kernel_settings& settings)
{
  #if 0
  double* sr = state_r();
  double* sw = state_w();
  for (int i_data = 0; i_data < n_dof*n_elem; ++i_data)
  {
    sw[i_data] = sr[i_data];
  }
  #endif
}

void Deformed_grid::execute_neighbor_derivative(int i_var, int i_dim, Kernel_settings& settings)
{
}

void Deformed_grid::execute_av_flux(Kernel_settings& settings)
{
}

void Deformed_grid::execute_local_av(int i_var, int i_dim, Kernel_settings& settings)
{
}

void Deformed_grid::execute_neighbor_av(int i_var, int i_dim, Kernel_settings& settings)
{
}

void Deformed_grid::connect(std::array<int, 2> i_elem, std::array<int, 2> i_dim,
                            std::array<bool, 2> is_positive)
{
  // combine vertices
  std::array<std::vector<int>, 2> vertex_inds;
  std::array<int, 2> strides;
  for (int i_side : {0, 1})
  {
    int stride = n_vertices/2;
    for (int i = 0; i < i_dim[i_side]; ++i) stride /= 2;
    strides[i_side] = stride;
    for (int i_vertex = 0; i_vertex < n_vertices; ++i_vertex)
    {
      if ((i_vertex/stride)%2 == int(is_positive[i_side]))
      {
        vertex_inds[i_side].push_back(i_vertex);
      }
    }
  }
  if ((is_positive[0] != is_positive[1]) && (i_dim[0] != i_dim[1]))
  {
    if (n_dim == 3)
    {
      int stride = strides[0];
      if (i_dim[0] < i_dim[1]) stride /= 2;
      for (int i : {0, 1}) std::swap(vertex_inds[1][i*2/stride], vertex_inds[1][i*2/stride + stride]);
    }
    else
    {
      std::swap(vertex_inds[1][0], vertex_inds[1][1]);
    }
  }
  if ((i_dim[0] == 0 && i_dim[1] == 2) || (i_dim[0] == 2 && i_dim[1] == 0))
  {
    std::swap(vertex_inds[1][1], vertex_inds[1][2]);
  }
  for (int i_vertex = 0; i_vertex < n_vertices/2; ++i_vertex)
  {
    Vertex& vert0 = deformed_element(i_elem[0]).vertex(vertex_inds[0][i_vertex]);
    Vertex& vert1 = deformed_element(i_elem[1]).vertex(vertex_inds[1][i_vertex]);
    vert0.eat(vert1);
  }

  // create connection object
  elem_cons.emplace_back();
  for (int i_side : {0, 1})
  {
    elem_cons.back().element[i_side] = elements[i_elem[i_side]].get();
    elem_cons.back().i_dim[i_side] = i_dim[i_side];
    elem_cons.back().is_positive[i_side] = is_positive[i_side];
  }
}

void Deformed_grid::purge_vertices()
{
  for (std::vector<Vertex::Non_transferable_ptr>::iterator current = vertices.begin();
       current < vertices.end(); ++current)
  {
    if (!*current) vertices.erase(current);
  }
}

void Deformed_grid::calc_vertex_relaxation()
{
  for (Vertex::Non_transferable_ptr& vertex : vertices)
  {
    if (vertex) vertex->calc_relax();
  }
}

void Deformed_grid::apply_vertex_relaxation()
{
  for (Vertex::Non_transferable_ptr& vertex : vertices)
  {
    if (vertex) vertex->apply_relax();
  }
}

std::vector<double> Deformed_grid::face_integral(Domain_func& integrand, int i_elem, int i_dim, bool is_positive)
{
  auto pos = get_pos(i_elem);
  auto row_weights = basis.node_weights();
  int stride = std::pow(basis.row_size, n_dim - 1 - i_dim);
  Deformed_element& elem = deformed_element(i_elem);
  double* stage = elem.stage(0);
  std::vector<double> total;
  for (int i_outer = 0; i_outer < n_qpoint/basis.row_size/stride; ++i_outer)
  {
    for (int i_inner = 0; i_inner < stride; ++i_inner)
    {
      int i_qpoint = (i_outer*basis.row_size + int(is_positive)*(basis.row_size - 1.))*stride + i_inner;
      std::vector<double> qpoint_pos;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        qpoint_pos.push_back(pos[i_qpoint + n_qpoint*j_dim]);
      }
      std::vector<double> qpoint_state;
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        qpoint_state.push_back(stage[i_qpoint + i_var*n_qpoint]);
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
          qpoint_jacobian(j_dim, k_dim) = elem.jacobian(j_dim, k_dim, i_qpoint);
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
        int stride_j = std::pow(basis.row_size, j_dim);
        weight *= row_weights((i_face/stride_j)%basis.row_size);
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
  for (Deformed_elem_wall wall : walls)
  {
    auto fi = face_integral(integrand, wall.i_elem, wall.i_dim, wall.is_positive);
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
  // FIXME: make this a kernel
  auto diff_mat = basis.diff_mat();
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> elem_pos = get_pos(i_elem);
    double* jac = elements[i_elem]->jacobian();
    for (int i_dim = 0, stride = n_qpoint/basis.row_size; i_dim < n_dim; ++i_dim, stride /= basis.row_size)
    {
      for (int i_outer = 0; i_outer < n_qpoint/(stride*basis.row_size); ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          for (int j_dim = 0; j_dim < n_dim; ++j_dim)
          {
            Eigen::VectorXd row_pos (basis.row_size);
            int row_start = i_outer*stride*basis.row_size + i_inner;
            for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint)
            {
              row_pos(i_qpoint) = elem_pos[j_dim*n_qpoint + row_start + i_qpoint*stride];
            }
            auto row_jacobian = diff_mat*row_pos;
            for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint)
            {
              jac[(j_dim*n_dim + i_dim)*n_qpoint + row_start + i_qpoint*stride] = row_jacobian(i_qpoint)/mesh_size;
            }
          }
        }
      }
    }
  }
}

std::string Deformed_grid::annotate(std::string file_name)
{
  return file_name + "_deformed";
}

}
