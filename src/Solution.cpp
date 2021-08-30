#include <limits>
#include <iostream>

#include <Solution.hpp>

namespace cartdg
{

#define FOR_ALL_GRIDS(code) for (Grid* grid : all_grids()) { code }

Solution::Solution(int n_var_arg, int n_dim_arg, int row_size_arg, double bms)
: n_var(n_var_arg), n_dim(n_dim_arg), base_mesh_size(bms), basis(row_size_arg) {}

Solution::~Solution() {}

Grid& Solution::get_grid(int order_added)
{
  if ((int)grids.size() <= order_added)
  {
    throw std::runtime_error("The requested Grid does not exist.");
  }
  else
  {
    return grids[order_added];
  }
}

void Solution::visualize(std::string file_prefix)
{
  char buffer [100];
  FOR_ALL_GRIDS
  (
    snprintf(buffer, 100, "%s_%.2e_%.2e", file_prefix.c_str(), grid->mesh_size, grid->time);
    grid->visualize(std::string(buffer));
  )
}

std::vector<double> Solution::integral()
{
  State_variables state_variables;
  return integral(state_variables);
}
std::vector<double> Solution::integral(Domain_func& integrand)
{
  if (grids.empty())
  {
    return std::vector<double> {};
  }
  else
  {
    std::vector<double> total;
    FOR_ALL_GRIDS
    (
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
    )
    return total;
  }
}

std::vector<double> Solution::surface_integral(Domain_func& integrand)
{
  if (def_grids.empty())
  {
    return std::vector<double> {};
  }
  else
  {
    std::vector<double> total;
    for (Deformed_grid grid : def_grids)
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
  for (Grid& g : grids)
  {
    all.push_back(&g);
  }
  for (Deformed_grid& g : def_grids)
  {
    all.push_back(&g);
  }
  return all;
}

double Solution::update(double cfl_by_stable_cfl)
{
  double dt = std::numeric_limits<double>::max();
  FOR_ALL_GRIDS // FIXME: incorporate_jacobian
  (
    dt = std::min<double>(dt, grid->stable_time_step(cfl_by_stable_cfl, kernel_settings));
  )
  for (int i_rk = 0; i_rk < 3; ++i_rk)
  {
    FOR_ALL_GRIDS
    (
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_local(kernel_settings);
    )
    FOR_ALL_GRIDS
    (
      kernel_settings.d_t_by_d_pos = dt/grid->mesh_size;
      grid->execute_neighbor(kernel_settings);
    )
    FOR_ALL_GRIDS
    (
      grid->execute_runge_kutta_stage();
    )
  }

  #if 1
  int n_iter = 4;
  double visc_dt = dt/n_iter;
  FOR_ALL_GRIDS // FIXME: incorporate_jacobian
  (
    kernel_settings.d_pos = grid->mesh_size;
    grid->execute_req_visc(kernel_settings);
  )
  FOR_ALL_GRIDS
  (
    kernel_settings.d_pos = grid->mesh_size;
    grid->execute_cont_visc(kernel_settings);
  )
  for (int i_iter = 0; i_iter < n_iter; ++i_iter)
  {
    for (int i_rk = 0; i_rk < 3; ++i_rk)
    {
      FOR_ALL_GRIDS
      (
        double* sr = grid->state_r();
        double* sw = grid->state_w();
        for (int i_data = 0; i_data < grid->n_elem*grid->n_dof; ++i_data)
        {
          sw[i_data] = sr[i_data];
        }
      )
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          FOR_ALL_GRIDS
          (
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_local_derivative(i_var, i_dim, kernel_settings);
          )
          FOR_ALL_GRIDS
          (
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_neighbor_derivative(i_var, i_dim, kernel_settings);
          )
          FOR_ALL_GRIDS
          (
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_av_flux(kernel_settings);
          )
          FOR_ALL_GRIDS
          (
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_local_av(i_var, i_dim, kernel_settings);
          )
          FOR_ALL_GRIDS
          (
            kernel_settings.d_t_by_d_pos = visc_dt/grid->mesh_size;
            kernel_settings.d_pos = grid->mesh_size;
            grid->execute_neighbor_av(i_var, i_dim, kernel_settings);
          )
        }
      }
      for (Grid& g : grids)
      {
        g.execute_runge_kutta_stage();
      }
    }
  }
  #endif
  FOR_ALL_GRIDS
  (
    grid->time += dt;
  )
  return dt;
}

void Solution::initialize(Spacetime_func& init_cond)
{
  FOR_ALL_GRIDS
  (
    double* state = grid->state_r();
    for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem)
    {
      std::vector<double> pos = grid->get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid->n_qpoint; ++i_qpoint)
      {
        std::vector<double> qpoint_pos;
        for (int i_dim = 0; i_dim < grid->n_dim; ++i_dim)
        {
          qpoint_pos.push_back(pos[i_qpoint + i_dim*grid->n_qpoint]);
        }
        auto qpoint_state = init_cond(qpoint_pos, grid->time);
        int qpoint_ind = i_elem*grid->n_dof + i_qpoint;
        for (int i_var = 0; i_var < grid->n_var; ++i_var)
        {
          state[qpoint_ind + i_var*grid->n_qpoint] = qpoint_state[i_var];
        }
      }
    }
  )
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
  grids.emplace_back(n_var, n_dim, n_elem, refined_mesh_size(ref_level), basis);
  Grid& g = grids.back();
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
  grids.emplace_back(n_var, n_dim, 0, refined_mesh_size(ref_level), basis);
}

void Solution::add_deformed_grid(int ref_level)
{
  def_grids.emplace_back(n_var, n_dim, 0, refined_mesh_size(ref_level), basis);
}

void Solution::auto_connect()
{
  for (Grid& grid : grids)
  {
    grid.auto_connect();
  }
}

void Solution::clear_neighbors()
{
  for (Grid& grid : grids)
  {
    grid.clear_neighbors();
  }
}

#undef FOR_ALL_GRIDS

}
