#include <limits>

#include <Solution.hpp>
#include <Tecplot_file.hpp>

namespace cartdg
{

Solution::Solution(int n_var_arg, int n_dim_arg, int row_size_arg, double bms)
: n_var(n_var_arg), n_dim(n_dim_arg), base_mesh_size(bms), basis(row_size_arg), time(0.)
{}

Solution::~Solution() {}

void Solution::visualize_field(std::string name)
{
  Tecplot_file file {name, n_dim, n_var, time};
  for (Grid* grid : all_grids())
  {
    if (grid->n_elem > 0)
    {
      grid->visualize_qpoints (file);
      grid->visualize_edges   (file);
      grid->visualize_interior(file);
    }
  }
}

void Solution::visualize_surface(std::string name)
{
  if (n_dim == 1) return; // 1D doesn't have visualizable surfaces
  Tecplot_file file {name, n_dim, n_var, time};
  for (Deformed_grid& grid : def_grids)
  {
    grid.visualize_surface(file);
  }
}

std::vector<double> Solution::integral()
{
  State_variables state_variables;
  return integral(state_variables);
}
std::vector<double> Solution::integral(Domain_func& integrand)
{
  auto grids = all_grids();
  if (grids.empty())
  {
    return std::vector<double> {};
  }
  else
  {
    std::vector<double> total;
    for(Grid* grid : all_grids())
    {
      auto grid_integral = grid->integral(integrand);
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
  for (Grid* grid : all_grids()) // FIXME: incorporate_jacobian
  {
    dt = std::min<double>(dt, grid->stable_time_step(cfl_by_stable_cfl, kernel_settings));
  }
  for (int i_rk = 0; i_rk < 3; ++i_rk)
  {
    for (Grid* grid : all_grids())
    {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_write_face(kernel_settings);
    }
    for (Grid* grid : all_grids())
    {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_neighbor(kernel_settings);
    }
    for (Grid* grid : all_grids())
    {
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_local(kernel_settings);
    }
    for (Deformed_grid& grid : def_grids) {
      grid.project_degenerate(kernel_settings.i_write);
    }
    for (Grid* grid : all_grids())
    {
      grid->execute_runge_kutta_stage();
    }
  }

  #if 0
  int n_iter = 4;
  double visc_dt = dt/n_iter;
  for (Grid* grid : all_grids()) // FIXME: incorporate_jacobian
  {
    kernel_settings.d_pos = grid->mesh_size;
    grid->execute_req_visc(kernel_settings);
  }
  for (Grid* grid : all_grids())
  {
    kernel_settings.d_pos = grid->mesh_size;
    grid->execute_cont_visc(kernel_settings);
  }
  for (int i_iter = 0; i_iter < n_iter; ++i_iter)
  {
    for (int i_rk = 0; i_rk < 3; ++i_rk)
    {
      for (Grid* grid : all_grids())
      {
        for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem)
        {
          Element& elem = grid->element(i_elem);
          double* stage_r = elem.stage(grid->i_stage_read());
          double* stage_w = elem.stage(grid->i_stage_write());
          for (int i_dof = 0; i_dof < grid->n_dof; ++i_dof)
          {
            stage_w[i_dof] = stage_r[i_dof];
          }
        }
      }
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          for (Grid* grid : all_grids())
          {
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_local_derivative(i_var, i_dim, kernel_settings);
          }
          for (Grid* grid : all_grids())
          {
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_neighbor_derivative(i_var, i_dim, kernel_settings);
          }
          for (Grid* grid : all_grids())
          {
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_av_flux(kernel_settings);
          }
          for (Grid* grid : all_grids())
          {
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_local_av(i_var, i_dim, kernel_settings);
          }
          for (Grid* grid : all_grids())
          {
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_neighbor_av(i_var, i_dim, kernel_settings);
          }
        }
      }
      for (Grid* grid : all_grids())
      {
        grid->execute_runge_kutta_stage();
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

void Solution::share_vertex_data(Element::shareable_value_access access_func)
{
  for (Grid* grid : all_grids()) {
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem) {
      grid->element(i_elem).push_shareable_value(access_func);
    }
  }
  for (Grid* grid : all_grids()) {
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem) {
      grid->element(i_elem).fetch_shareable_value(access_func);
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
