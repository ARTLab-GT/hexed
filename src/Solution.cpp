#include <limits>

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
    throw "The requested Grid does not exist.";
  }
  else
  {
    return grids[order_added];
  }
}

void Solution::visualize(std::string file_prefix)
{
  char buffer [100];
  for (Grid* grid : all_grids)
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
    for (Grid& grid : grids)
    {
      auto grid_integral = grid.integral(integrand);
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

double Solution::update(double cfl_by_stable_cfl)
{
  Local_kernel local = get_local_kernel();
  Neighbor_kernel neighbor = get_neighbor_kernel();
  Max_char_speed_kernel max_char_speed = get_max_char_speed_kernel();
  Fbc_kernel fbc = get_fbc_kernel();
  double dt = std::numeric_limits<double>::max();
  for (Grid& g : grids)
  {
    double cfl = cfl_by_stable_cfl*g.get_stable_cfl();
    dt = std::min<double>(dt, cfl*g.mesh_size/max_char_speed(g.state_r(), g.n_elem,
                                                             kernel_settings));
  }
  for (Grid& g : grids)
  {
    do
    {
      kernel_settings.d_t_by_d_pos = dt/g.mesh_size;
      local(g.state_r(), g.state_w(), g.n_elem, g.basis.diff_mat(), kernel_settings);
      neighbor(g.neighbor_connections_r().data(), g.neighbor_connections_w().data(), 
               g.n_neighb_con().data(), g.basis.node_weights(), kernel_settings);
      {
        auto weights = g.basis.node_weights();
        fbc(g.fit_bound_conds, g.state_r(), g.state_w(), weights(0), kernel_settings);
      }
    }
    while (!g.execute_runge_kutta_stage());
    g.time += dt;
  }
  return dt;
}

void Solution::initialize(Spacetime_func& init_cond)
{
  for (Grid& g : grids)
  {
    double* state = g.state_r();
    for (int i_elem = 0; i_elem < g.n_elem; ++i_elem)
    {
      std::vector<double> pos = g.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < g.n_qpoint; ++i_qpoint)
      {
        std::vector<double> qpoint_pos;
        for (int i_dim = 0; i_dim < g.n_dim; ++i_dim)
        {
          qpoint_pos.push_back(pos[i_qpoint + i_dim*g.n_qpoint]);
        }
        auto qpoint_state = init_cond(qpoint_pos, g.time);
        int qpoint_ind = i_elem*g.n_dof + i_qpoint;
        for (int i_var = 0; i_var < g.n_var; ++i_var)
        {
          state[qpoint_ind + i_var*g.n_qpoint] = qpoint_state[i_var];
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
  all_grids.push_back(&g);
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
  all_grids.push_back(&grids.back());
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
