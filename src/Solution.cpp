#include <limits>
#include <iostream>

#include <Solution.hpp>

namespace cartdg
{

Solution::Solution(int n_var_arg, int n_dim_arg, int rank_arg, double bms)
: n_var(n_var_arg), n_dim(n_dim_arg), base_mesh_size(bms), basis(rank_arg) {}

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
  for (Grid* grid : all_grids())
  {
    snprintf(buffer, 100, "%s_%.2e_%.2e", file_prefix.c_str(), grid->mesh_size, grid->time);
    grid->visualize(std::string(buffer));
  }
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
    for (Grid* grid : all_grids())
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
  auto local = get_local_kernel();
  auto local_deformed = get_local_deformed_kernel();
  auto neighbor = get_neighbor_kernel();
  auto neighbor_deformed = get_neighbor_deformed_kernel();
  auto nonpen = get_nonpen_kernel();
  auto max_char_speed = get_max_char_speed_kernel();
  auto physical_step = get_physical_step_kernel();
  auto restrict_step = get_restrict_step_kernel();
  auto fbc = get_gbc_kernel();
  double dt = std::numeric_limits<double>::max();
  for (Grid* g : all_grids()) // FIXME: incorporate jacobian
  {
    double cfl = cfl_by_stable_cfl*g->get_stable_cfl();
    dt = std::min<double>(dt, cfl*g->mesh_size/max_char_speed(g->state_r(), g->n_elem,
                                                              kernel_settings));
  }
  for (int i_rk = 0; i_rk < 3; ++i_rk)
  {
    double step = 1.;
    for (Grid& g : grids)
    {
      kernel_settings.d_t_by_d_pos = dt/g.mesh_size;
      local(g.state_r(), g.state_w(), g.n_elem, g.basis, kernel_settings);
      neighbor(g.neighbor_connections_r().data(), g.neighbor_connections_w().data(), 
               g.n_neighb_con().data(), g.basis.node_weights(), kernel_settings);
      {
        auto weights = g.basis.node_weights();
        fbc(g.ghost_bound_conds, g.state_r(), g.state_w(), weights(0), kernel_settings);
      }
      double restricted_step = physical_step(g.state_r(), g.state_w(), g.n_elem, kernel_settings);
      if (restricted_step < step)
      {
        step = restricted_step;
        printf("\nTime step restricted to prevent nonphysical thermodynamic quantities!\n");
      }
    }
    for (Deformed_grid& g : def_grids)
    {
      kernel_settings.d_t_by_d_pos = dt/g.mesh_size;
      local_deformed(g.state_r(), g.state_w(), g.jacobian.data(), g.n_elem,
                     g.basis, kernel_settings);
      neighbor_deformed(g.state_connections_r(), g.state_connections_w(), g.jacobian_neighbors.data(),
                        g.neighbor_axes.data(), g.neighbor_is_positive.data(), g.neighbor_storage[0].size()/2, g.basis.node_weights(), kernel_settings);
      {
        auto weights = g.basis.node_weights();
        fbc(g.ghost_bound_conds, g.state_r(), g.state_w(), weights(0), kernel_settings);
      }
      nonpen(g.state_r(), g.state_w(), g.jacobian.data(), g.i_elem_wall.data(), g.i_dim_wall.data(), g.is_positive_wall.data(), g.i_elem_wall.size(), g.basis.node_weights()(0), kernel_settings);
      double restricted_step = physical_step(g.state_r(), g.state_w(), g.n_elem, kernel_settings);
      if (restricted_step < step)
      {
        step = restricted_step;
        printf("\nTime step restricted to prevent nonphysical thermodynamic quantities!\n");
      }
    }
    for (Grid* g : all_grids())
    {
      restrict_step(g->state_r(), g->state_w(), g->n_elem, step, kernel_settings);
      g->execute_runge_kutta_stage();
    }
  }
  int n_iter = 30;
  for (int i_iter = 0; i_iter < n_iter; ++i_iter)
  {
    dt = std::numeric_limits<double>::max();
    for (Grid* g : all_grids()) // FIXME: incorporate jacobian
    {
      double cfl = cfl_by_stable_cfl*g->get_stable_cfl();
      dt = std::min<double>(dt, cfl*g->mesh_size/max_char_speed(g->state_r(), g->n_elem,
                                                                kernel_settings));
    }
    for (int i_rk = 0; i_rk < 3; ++i_rk)
    {
      for (Grid& g : grids)
      {
        kernel_settings.d_t_by_d_pos = dt/g.mesh_size;
        double* sr = g.state_r();
        double* sw = g.state_w();
        for (int i_data = 0; i_data < g.n_dof*g.n_elem; ++i_data)
        {
          sw[i_data] = sr[i_data];
        }
        for (int i_axis = 0; i_axis < n_dim; ++i_axis)
        {
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            get_derivative_kernel()(sr, g.derivs.data(), g.n_elem, i_var, i_axis, basis, kernel_settings);
            int n_con = g.n_neighb_con()[i_axis];
            get_jump_kernel()(g.neighbor_connections_r()[i_axis],
                              g.deriv_neighbor_connections()[i_axis], n_con, i_var, i_axis, basis.node_weights(), kernel_settings);
            for (int i_data = 0; i_data < g.n_qpoint*g.n_elem; ++i_data)
            {
              g.derivs[i_data] *= 1.e-4;
            }
            get_viscous_local_kernel()(g.derivs.data(), sw, g.n_elem, i_var, i_axis, basis, kernel_settings);
            get_viscous_neighbor_kernel()(g.deriv_neighbor_connections()[i_axis],
                                          g.neighbor_connections_w()[i_axis], n_con, i_var, i_axis, basis.node_weights(), kernel_settings);
            get_jump_gbc_kernel()(g.ghost_bound_conds, g.derivs.data(), g.state_w(), i_var, i_axis, basis.node_weights()(0), kernel_settings);
          }
        }
      }
      for (Grid& g : grids)
      {
        g.execute_runge_kutta_stage();
      }
    }
  }
  for (Grid* g : all_grids())
  {
    g->time += dt;
  }
  return dt;
}

void Solution::initialize(Spacetime_func& init_cond)
{
  for (Grid* g : all_grids())
  {
    double* state = g->state_r();
    for (int i_elem = 0; i_elem < g->n_elem; ++i_elem)
    {
      std::vector<double> pos = g->get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < g->n_qpoint; ++i_qpoint)
      {
        std::vector<double> qpoint_pos;
        for (int i_dim = 0; i_dim < g->n_dim; ++i_dim)
        {
          qpoint_pos.push_back(pos[i_qpoint + i_dim*g->n_qpoint]);
        }
        auto qpoint_state = init_cond(qpoint_pos, g->time);
        int qpoint_ind = i_elem*g->n_dof + i_qpoint;
        for (int i_var = 0; i_var < g->n_var; ++i_var)
        {
          state[qpoint_ind + i_var*g->n_qpoint] = qpoint_state[i_var];
        }
      }
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

}
