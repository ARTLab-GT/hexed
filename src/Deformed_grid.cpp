#include <Deformed_grid.hpp>
#include <Regular_grid.hpp>
#include <Tecplot_file.hpp>
#include <Physical_basis.hpp>
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
  Deformed_elem_wall wall {{elements[i_elem].get(), i_dim, is_positive_face}, i_elem};
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

void Deformed_grid::project_degenerate(int i_stage)
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    Deformed_element& elem {deformed_element(i_elem)};
    if (elem.degenerate) {
      Physical_basis phys {n_dim, basis.row_size, get_pos(i_elem)};
      Eigen::VectorXd weights (n_qpoint);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        double product = 1.;
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          int i_row = (i_qpoint/custom_math::pow(basis.row_size, n_dim - 1 - i_dim))%basis.row_size;
          product *= basis.node_weights()(i_row);
        }
        product *= elem.jacobian_determinant(i_qpoint);
        weights[i_qpoint] = product;
      }
      Eigen::Map<Eigen::MatrixXd> state (elem.stage(i_stage), n_qpoint, n_var);
      state = phys.projection(state, weights);
    }
  }
}

void Deformed_grid::connect(std::array<int, 2> i_elem, std::array<int, 2> i_dim,
                            std::array<bool, 2> is_positive)
{
  // create connection object
  std::array<Face_index, 2> face_inds;
  for (int i_side : {0, 1})
  {
    face_inds[i_side].element = elements[i_elem[i_side]].get();
    face_inds[i_side].i_dim = i_dim[i_side];
    face_inds[i_side].is_positive = is_positive[i_side];
  }
  elem_cons.emplace_back(face_inds);
  Deformed_elem_con& con = elem_cons.back();

  // get vertices involved
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

  // reorder as necessary
  if (con.flip_tangential()) {
    for (int i_vert = 0; i_vert < n_vertices/2; ++i_vert) {
      // this arithmetic is equivalent to swapping vertices along dimension i_dim[0]
      int coord = (vertex_inds[1][i_vert]/strides[0])%2;
      vertex_inds[1][i_vert] += ((coord + 1)%2 - coord)*strides[0];
    }
  }
  if (con.transpose()) std::swap(vertex_inds[1][1], vertex_inds[1][2]);

  // combine vertices
  for (int i_vertex = 0; i_vertex < n_vertices/2; ++i_vertex)
  {
    Vertex& vert0 = deformed_element(i_elem[0]).vertex(vertex_inds[0][i_vertex]);
    Vertex& vert1 = deformed_element(i_elem[1]).vertex(vertex_inds[1][i_vertex]);
    vert0.eat(vert1);
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

std::vector<double> Deformed_grid::face_integral(Surface_func& integrand, int i_elem, int i_dim, bool is_positive)
{
  auto pos = get_pos(i_elem);
  auto row_weights = basis.node_weights();
  Deformed_element& elem = deformed_element(i_elem);
  double* stage = elem.stage(0);
  double* jacobian {elem.jacobian()};
  int nfq {n_qpoint/basis.row_size};
  Eigen::MatrixXd boundary {basis.boundary().row(is_positive)};
  Eigen::MatrixXd face_pos {nfq, n_dim};
  Eigen::MatrixXd face_state {nfq, n_var};
  Eigen::MatrixXd face_jac {nfq, n_dim*n_dim};
  for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
    Eigen::Map<Eigen::VectorXd> pos_vec {pos.data() + j_dim*n_qpoint, n_qpoint};
    face_pos.col(j_dim) = custom_math::dimension_matvec(boundary, pos_vec, i_dim);
  }
  for (int i_var = 0; i_var < n_var; ++i_var) {
    Eigen::Map<Eigen::VectorXd> state_vec {stage + i_var*n_qpoint, n_qpoint};
    face_state.col(i_var) = custom_math::dimension_matvec(boundary, state_vec, i_dim);
  }
  for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) {
    Eigen::Map<Eigen::VectorXd> jac_vec {jacobian + i_jac*n_qpoint, n_qpoint};
    face_jac.col(i_jac) = custom_math::dimension_matvec(boundary, jac_vec, i_dim);
  }
  std::vector<double> total;

  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint)
  {
    std::vector<double> qpoint_pos;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      qpoint_pos.push_back(face_pos(i_qpoint, j_dim));
    }
    std::vector<double> qpoint_state;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      qpoint_state.push_back(face_state(i_qpoint, i_var));
    }
    Eigen::MatrixXd qpoint_jacobian (n_dim, n_dim);
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
        qpoint_jacobian(j_dim, k_dim) = face_jac(i_qpoint, j_dim*n_dim + k_dim);
      }
      qpoint_jacobian(j_dim, i_dim) = 0.;
    }
    std::vector<double> normal;
    double face_jac_det {0.};
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      qpoint_jacobian(j_dim, i_dim) = 1.;
      double normal_j {qpoint_jacobian.determinant()*(1 - 2*is_positive)};
      qpoint_jacobian(j_dim, i_dim) = 0.;
      normal.push_back(normal_j);
      face_jac_det += normal_j*normal_j;
    }
    face_jac_det = std::sqrt(face_jac_det);
    auto qpoint_integrand = integrand(qpoint_pos, time, qpoint_state, normal);
    if (total.size() < qpoint_integrand.size()) {
      total.resize(qpoint_integrand.size());
    }
    double weight = 1.;
    for (int j_dim = 0; j_dim < n_dim - 1; ++j_dim) {
      int stride_j = custom_math::pow(basis.row_size, j_dim);
      weight *= row_weights((i_qpoint/stride_j)%basis.row_size);
    }
    for (int i_var = 0; i_var < int(total.size()); ++i_var) {
      total[i_var] += weight*face_jac_det*qpoint_integrand[i_var];
    }
  }
  for (int i_var = 0; i_var < int(total.size()); ++i_var) {
    total[i_var] *= custom_math::pow(mesh_size, n_dim - 1);
  }
  return total;
}

std::vector<double> Deformed_grid::surface_integral(Surface_func& integrand)
{
  std::vector<double> total;
  for (Deformed_elem_wall wall : walls)
  {
    auto face_ind = wall.face_index();
    auto fi = face_integral(integrand, wall.i_elem(), face_ind.i_dim, face_ind.is_positive);
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

void Deformed_grid::visualize_connections(std::string file_name)
{
  int n_con = elem_cons.size();
  Tecplot_file file {file_name, n_dim, 0, time};
  Tecplot_file::Line_segments segs {file, n_con, 2, "connections"};
  for (int i_con = 0; i_con < n_con; ++i_con)
  {
    Deformed_elem_con& con {elem_cons[i_con]};
    std::vector<double> centers (2*n_dim, 0.);
    for (int i_side : {0, 1}) {
      for (int i_vert = 0; i_vert < n_vertices; ++i_vert) {
        auto pos {con.face_index(i_side).element->vertex(i_vert).pos};
        double weight {0.};
        if (i_vert == 0) weight += i_con;
        else if (i_vert == n_vertices - 1) weight += n_con - i_con;
        weight = (weight/n_con + 1.)/(n_vertices + 1.);
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          centers[2*i_dim + i_side] += pos[i_dim]*weight;
        }
      }
    }
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      double avg {(centers[2*i_dim] + centers[2*i_dim + 1])/2.};
      for (int i_side : {0, 1}) {
        centers[2*i_dim + i_side] = (centers[2*i_dim + i_side] + avg)/2.;
      }
    }
    segs.write(centers.data(), nullptr);
  }
}

void Deformed_grid::visualize_surface(Tecplot_file& file, int n_sample)
{
  if (n_dim == 1) throw std::runtime_error("cannot visualize surfaces in 1D");
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  Eigen::MatrixXd boundary {basis.boundary()};
  const int n_block {custom_math::pow(n_sample, n_dim - 1)};
  for (Deformed_elem_wall wall : walls)
  {
    auto fi = wall.face_index();
    std::vector<double> pos = get_pos(wall.i_elem());
    Eigen::VectorXd interp_pos {n_block*n_dim};
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      Eigen::Map<Eigen::VectorXd> qpoint_pos (pos.data() + i_dim*n_qpoint, n_qpoint);
      auto face {custom_math::dimension_matvec(boundary.row(fi.is_positive), qpoint_pos, fi.i_dim)};
      interp_pos.segment(i_dim*n_block, n_block) = custom_math::hypercube_matvec(interp, face);
    }
    double* state = element(wall.i_elem()).stage(0);
    Eigen::VectorXd interp_state {n_block*n_var};
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      Eigen::Map<Eigen::VectorXd> var (state + i_var*n_qpoint, n_qpoint);
      auto face {custom_math::dimension_matvec(boundary.row(fi.is_positive), var, fi.i_dim)};
      interp_state.segment(i_var*n_block, n_block) = custom_math::hypercube_matvec(interp, face);
    }
    Tecplot_file::Structured_block zone {file, n_sample, "face_interior", n_dim - 1};
    zone.write(interp_pos.data(), interp_state.data());
  }
}

void Deformed_grid::calc_jacobian()
{
  // compute jacobian of each element
  auto diff_mat = basis.diff_mat();
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    std::vector<double> elem_pos = get_pos(i_elem);
    double* jac = elements[i_elem]->jacobian();
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      Eigen::Map<Eigen::VectorXd> dim_pos (elem_pos.data() + i_dim*n_qpoint, n_qpoint);
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        Eigen::Map<Eigen::VectorXd> jac_entry (jac + (i_dim*n_dim + j_dim)*n_qpoint, n_qpoint);
        jac_entry = custom_math::dimension_matvec(diff_mat, dim_pos, j_dim)/mesh_size;
      }
    }
  }
}

std::string Deformed_grid::annotate(std::string file_name)
{
  return file_name + "_deformed";
}

}
