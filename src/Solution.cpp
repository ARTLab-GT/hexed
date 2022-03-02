#include <limits>

#include <Solution.hpp>
#include <Tecplot_file.hpp>
#include <math.hpp>

namespace cartdg
{

Solution::Solution(int n_var_arg, int n_dim_arg, int row_size_arg, double bms)
: n_var(n_var_arg), n_dim(n_dim_arg), base_mesh_size(bms), basis(row_size_arg), time(0.)
{}

Solution::~Solution() {}

void Solution::visualize_field(Qpoint_func& func, std::string name)
{
  const int n_vis = func.n_var(n_dim); // number of variables to visualize
  std::vector<std::string> var_names;
  for (int i_vis = 0; i_vis < n_vis; ++i_vis) var_names.push_back(func.variable_name(i_vis));
  Tecplot_file file {name, n_dim, var_names, time};
  const int n_sample = 20;
  for (Grid* grid : all_grids())
  {
    if (grid->n_elem > 0)
    {
      const int n_elem = grid->n_elem;
      const int n_qpoint = grid->n_qpoint;
      const int n_corners {custom_math::pow(2, n_dim - 1)};
      const int nfqpoint = n_qpoint/basis.row_size;
      const int n_block {custom_math::pow(n_sample, n_dim)};
      Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};

      for (int i_elem = 0; i_elem < n_elem; ++i_elem)
      {
        std::vector<double> pos = grid->get_pos(i_elem);
        std::vector<double> to_vis (n_qpoint*n_vis);
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          auto qpoint_vis = func(*grid, i_elem, i_qpoint);
          for (int i_vis = 0; i_vis < n_vis; ++i_vis) {
            to_vis[i_vis*n_qpoint + i_qpoint] = qpoint_vis[i_vis];
          }
        }
        // note: each visualization stage is enclosed in `{}` to ensure that only one
        // `Tecplot_file::Zone` is alive at a time

        // visualize edges
        if (n_dim > 1) // 1D elements don't really have edges
        {
          Tecplot_file::Line_segments edges {file, n_dim*n_corners, n_sample, "edges"};
          Eigen::MatrixXd boundary {basis.boundary()};
          for (int i_dim = 0; i_dim < n_dim; ++i_dim)
          {
            const int stride {custom_math::pow(basis.row_size, n_dim - 1 - i_dim)};
            const int n_outer {n_qpoint/stride/basis.row_size};

            auto extract_edge = [=](double* data, int n)
            {
              Eigen::MatrixXd edge {n_sample, n_corners*n};
              for (int i = 0; i < n; ++i) {
                Eigen::MatrixXd edge_qpoints {basis.row_size, n_corners};
                for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint) {
                  Eigen::VectorXd qpoint_slab {nfqpoint};
                  for (int i_outer = 0; i_outer < n_outer; ++i_outer) {
                    for (int i_inner = 0; i_inner < stride; ++i_inner) {
                      qpoint_slab[i_outer*stride + i_inner] = data[i*n_qpoint + i_qpoint*stride + i_outer*stride*basis.row_size + i_inner];
                    }
                  }
                  edge_qpoints.row(i_qpoint) = custom_math::hypercube_matvec(boundary, qpoint_slab);
                }
                for (int i_corner = 0; i_corner < n_corners; ++i_corner) {
                  edge.col(i_corner*n + i) = interp*edge_qpoints.col(i_corner);
                }
              }
              return edge;
            };

            Eigen::MatrixXd edge_pos {extract_edge(pos.data(), n_dim)};
            Eigen::MatrixXd edge_state {extract_edge(to_vis.data(), n_vis)};
            for (int i_corner = 0; i_corner < n_corners; ++i_corner) {
              edges.write(edge_pos.data() + i_corner*n_dim*n_sample, edge_state.data() + i_corner*n_vis*n_sample);
            }
          }
        }

        { // visualize quadrature points
          Tecplot_file::Structured_block qpoints {file, basis.row_size, "element_qpoints"};
          qpoints.write(pos.data(), to_vis.data());
        }

        { // visualize interior (that is, quadrature point data interpolated to a fine mesh of sample points)
          Tecplot_file::Structured_block interior {file, n_sample, "element_interior"};
          Eigen::VectorXd interp_pos {n_block*n_dim};
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            Eigen::Map<Eigen::VectorXd> qpoint_pos (pos.data() + i_dim*n_qpoint, n_qpoint);
            interp_pos.segment(i_dim*n_block, n_block) = custom_math::hypercube_matvec(interp, qpoint_pos);
          }
          Eigen::VectorXd interp_state {n_block*n_vis};
          for (int i_var = 0; i_var < n_vis; ++i_var) {
            Eigen::Map<Eigen::VectorXd> var (to_vis.data() + i_var*n_qpoint, n_qpoint);
            interp_state.segment(i_var*n_block, n_block) = custom_math::hypercube_matvec(interp, var);
          }
          interior.write(interp_pos.data(), interp_state.data());
        }
      }
    }
  }
}

void Solution::visualize_surface(std::string name)
{
  if (n_dim == 1) return; // 1D doesn't have visualizable surfaces
  int n_vis = n_var; // FIXME
  std::vector<std::string> var_names;
  for (int i_vis = 0; i_vis < n_vis; ++i_vis) var_names.push_back("state" + std::to_string(i_vis));
  Tecplot_file file {name, n_dim, var_names, time};
  for (Deformed_grid& grid : def_grids) {
    grid.visualize_surface(file);
  }
}

std::vector<double> Solution::integral(Qpoint_func& integrand)
{
  auto grids = all_grids();
  if (grids.empty())
  {
    return std::vector<double> {};
  }
  else
  {
    // compute `n_dim`-dimensional quadrature weights
    int n_qpoint = custom_math::pow(basis.row_size, n_dim);
    Eigen::VectorXd weights (n_qpoint);
    Eigen::VectorXd weights_1d = basis.node_weights();
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) weights(i_qpoint) = 1.;
    for (int stride = n_qpoint/basis.row_size, n_rows = 1; n_rows < n_qpoint;
         stride /= basis.row_size, n_rows *= basis.row_size)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer) {
        for (int i_inner = 0; i_inner < stride; ++i_inner) {
          for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint) {
            weights((i_outer*basis.row_size + i_qpoint)*stride + i_inner) *= weights_1d(i_qpoint);
          }
        }
      }
    }

    // compute integral
    std::vector<double> total;
    for(Grid* grid : all_grids()) {
      double mesh_size_factor = custom_math::pow(grid->mesh_size, n_dim);
      for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem) {
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          std::vector<double> point_integrand = integrand(*grid, i_elem, i_qpoint);
          if (total.size() < point_integrand.size()) total.resize(point_integrand.size(), 0.);
          double jac_det = grid->element(i_elem).jacobian_determinant(i_qpoint);
          for (int i_var = 0; i_var < int(point_integrand.size()); ++i_var) {
            total[i_var] += point_integrand[i_var]*weights[i_qpoint]*jac_det*mesh_size_factor;
          }
        }
      }
    }
    return total;
  }
}

std::vector<double> Solution::surface_integral(Surface_func& integrand)
{
  if (def_grids.empty())
  {
    return std::vector<double> {};
  }
  else
  {
    std::vector<double> total;
    for (Deformed_grid& grid : def_grids)
    {
      auto grid_integral = grid.surface_integral(integrand);
      int size = grid_integral.size();
      if (int(total.size()) < size)
      {
        total.resize(size);
      }
      for (int i_var = 0; i_var < size; ++i_var)
      {
        total[i_var] += grid_integral[i_var];
      }
    }
    return total;
  }
}

std::vector<Grid*> Solution::all_grids()
{
  std::vector<Grid*> all;
  for (Grid& g : reg_grids)
  {
    all.push_back(&g);
  }
  for (Grid& g : def_grids)
  {
    all.push_back(&g);
  }
  return all;
}

double Solution::update(double cfl_by_stable_cfl)
{
  double dt = std::numeric_limits<double>::max();
  for (Grid* grid : all_grids()) {
    dt = std::min<double>(dt, grid->stable_time_step(cfl_by_stable_cfl, kernel_settings));
  }
  for (int i_rk = 0; i_rk < 3; ++i_rk)
  {
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_write_face(kernel_settings);
    }
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_neighbor(kernel_settings);
    }
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_local(kernel_settings);
    }
    for (Deformed_grid& grid : def_grids) {
      grid.project_degenerate(kernel_settings.i_write);
    }
    for (Grid* grid : all_grids()) {
      grid->execute_runge_kutta_stage();
    }
  }
  #if 1
  for (int i = 0; i < 10; ++i) {
  for (int i_var = 0; i_var < n_var; ++i_var)
  {
    for (Grid* grid : all_grids()) {
      kernel_settings.d_pos = grid->mesh_size;
      grid->execute_req_visc(kernel_settings);
    }
    share_vertex_data(&Element::viscosity);
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = 1./grid->mesh_size; // for gradient `d_t_by_d_pos` is interpeted as mesh spacing
      grid->execute_write_face_gradient(i_var, kernel_settings);
    }
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = 1./grid->mesh_size;
      grid->execute_neighbor_gradient(i_var, kernel_settings);
    }
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = 1./grid->mesh_size;
      grid->execute_local_gradient(i_var, kernel_settings);
    }
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_write_face_av(i_var, kernel_settings);
    }
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_neighbor_av(i_var, kernel_settings);
    }
    for (Grid* grid : all_grids()) {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_local_av(i_var, kernel_settings);
    }
  }
  for (Deformed_grid& grid : def_grids) {
    grid.project_degenerate(kernel_settings.i_read);
  }
  }
  #endif
  #if 0
  for (Grid* grid : all_grids()) {
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem) {
      double* stage_read  = grid->element(i_elem).stage(kernel_settings.i_read);
      double* stage_write = grid->element(i_elem).stage(kernel_settings.i_write);
      for (int i_dof = 0; i_dof < grid->n_var*grid->n_qpoint; ++i_dof) {
        stage_write[i_dof] = stage_read[i_dof];
      }
    }
  }
  #endif

  time += dt;
  for (Grid* grid : all_grids())
  {
    grid->time = time;
  }
  return dt;
}

void Solution::initialize(Spacetime_func& init_cond)
{
  for (Grid* grid : all_grids())
  {
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem)
    {
      std::vector<double> pos = grid->get_pos(i_elem);
      double* elem_state = grid->element(i_elem).stage(0);
      for (int i_qpoint = 0; i_qpoint < grid->n_qpoint; ++i_qpoint)
      {
        std::vector<double> qpoint_pos;
        for (int i_dim = 0; i_dim < grid->n_dim; ++i_dim) {
          qpoint_pos.push_back(pos[i_qpoint + i_dim*grid->n_qpoint]);
        }
        auto qpoint_state = init_cond(qpoint_pos, grid->time);
        for (int i_var = 0; i_var < grid->n_var; ++i_var) {
          elem_state[i_var*grid->n_qpoint + i_qpoint] = qpoint_state[i_var];
        }
      }
    }
  }
}

void Solution::share_vertex_data(Element::shareable_value_access access_func, Vertex::reduction reduce)
{
  for (Grid* grid : all_grids()) {
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem) {
      grid->element(i_elem).push_shareable_value(access_func);
    }
  }
  for (Grid* grid : all_grids()) {
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem) {
      grid->element(i_elem).fetch_shareable_value(access_func, reduce);
    }
  }
  for (Grid* grid : all_grids()) {
    grid->match_hanging(access_func);
  }
}

void Solution::set_local_time_step()
{
  // set time step to be continuous at vertices
  share_vertex_data(&Element::vertex_time_step_scale, Vertex::vector_min);
  // interpolate to quadrature points
  Eigen::MatrixXd lin_interp {basis.row_size, 2};
  for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint) {
    double node {basis.node(i_qpoint)};
    lin_interp(i_qpoint, 0) = 1. - node;
    lin_interp(i_qpoint, 1) = node;
  }
  for (Grid* grid : all_grids()) {
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem) {
      Element& elem {grid->element(i_elem)};
      Eigen::Map<Eigen::VectorXd> vert_tss (elem.vertex_time_step_scale(), grid->n_vertices);
      Eigen::Map<Eigen::VectorXd> qpoint_tss (elem.time_step_scale(), grid->n_qpoint);
      qpoint_tss = custom_math::hypercube_matvec(lin_interp, vert_tss);
    }
  }
}

double Solution::refined_mesh_size(int ref_level)
{
  double mesh_size = base_mesh_size;
  for (int i_rl = 0; i_rl < ref_level; ++i_rl) mesh_size /= 2;
  return mesh_size;
}

void Solution::add_block_grid(int ref_level, std::vector<int> lower_corner,
                                             std::vector<int> upper_corner)
{
  int n_elem = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_elem *= upper_corner[i_dim] - lower_corner[i_dim];
  }
  reg_grids.emplace_back(n_var, n_dim, n_elem, refined_mesh_size(ref_level), basis);
  Regular_grid& g = reg_grids.back();
  for (int i_dim = 0, stride = 1; i_dim < n_dim; ++i_dim)
  {
    int row_size = upper_corner[i_dim] - lower_corner[i_dim];
    for (int i_outer = 0; i_outer < n_elem/(stride*row_size); ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        for (int i_row = 0; i_row < row_size; ++i_row)
        {
          g.pos[n_dim*(i_row*stride + i_outer*stride*row_size + i_inner) + i_dim]
          = i_row + lower_corner[i_dim];
        }
      }
    }
    stride *= row_size;
  }
}

void Solution::add_block_grid(int ref_level)
{
  int n_div = 1;
  for (int i = 0; i < ref_level; ++i) n_div *= 2;
  std::vector<int> lc; std::vector<int> uc;
  for (int i = 0; i < n_dim; ++i)
  {
    lc.push_back(0); uc.push_back(n_div);
  }
  add_block_grid(ref_level, lc, uc);
}

void Solution::add_empty_grid(int ref_level)
{
  reg_grids.emplace_back(n_var, n_dim, 0, refined_mesh_size(ref_level), basis);
}

void Solution::add_deformed_grid(int ref_level)
{
  def_grids.emplace_back(n_var, n_dim, 0, refined_mesh_size(ref_level), basis);
}

void Solution::auto_connect()
{
  for (Regular_grid& grid : reg_grids)
  {
    grid.auto_connect();
  }
}

#undef FOR_ALL_GRIDS

}
