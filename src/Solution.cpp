#include <Solution.hpp>
#include <Initializer.hpp>

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

void Solution::update(double dt)
{
  Local_kernel local = get_local_kernel();
  Neighbor_kernel neighbor = get_neighbor_kernel();
  for (Grid& g : grids)
  {
    // FIXME: replace this with a CFL-based time step
    double d_t_by_d_x = dt/g.mesh_size;
    local(g.state_r(), g.state_w(), g.n_elem, g.basis.diff_mat(), d_t_by_d_x, 1.4);
    neighbor(g.neighbor_connections_r().data(), g.neighbor_connections_w().data(), 
             g.n_neighb_con().data(), g.basis.node_weights(), d_t_by_d_x, 1.4);
    g.time += dt;
    g.iter++;
  }
}

void Solution::update()
{
  update(1.);
}

void Solution::initialize(Initializer& init)
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
        std::vector<double> mmtm = init.momentum(qpoint_pos);
        std::vector<double> scalar = init.scalar_state(qpoint_pos);
        int qpoint_ind = i_elem*g.n_dof + i_qpoint;
        for (int i_dim = 0; i_dim < g.n_dim; ++i_dim)
        {
          state[qpoint_ind + i_dim*g.n_qpoint] = mmtm[i_dim];
        }
        for (int i_var = g.n_dim; i_var < g.n_var; ++i_var)
        {
          state[qpoint_ind + i_var*g.n_qpoint] = scalar[i_var - g.n_dim];
        }
      }
    }
  }
}

void Solution::add_block_grid(int ref_level, std::vector<int> lower_corner,
                                             std::vector<int> upper_corner)
{
  int n_elem = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_elem *= upper_corner[i_dim] - lower_corner[i_dim];
  }
  double mesh_size = base_mesh_size;
  for (int i_rl = 0; i_rl < ref_level; ++i_rl) mesh_size /= 2;
  grids.push_back(Grid (n_var, n_dim, n_elem, mesh_size, basis));
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

void Solution::auto_connect()
{
  for (Grid& grid : grids)
  {
    grid.auto_connect();
  }
}
