#include "Solution.hpp"

Solution::Solution(int n_var_arg, int n_dim_arg, int rank_arg, double bms)
: n_var(n_var_arg), n_dim(n_dim_arg), basis(rank_arg), base_mesh_size(bms) {}

Solution::~Solution() {}

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
  
  Grid g (n_var, n_dim, n_elem, mesh_size, basis);
  grids.push_back(g);
}

Grid& Solution::get_grid(int order_added)
{
  if (grids.size() <= order_added)
  {
    throw "The requested grid does not exist.";
  }
  else
  {
    return grids[order_added];
  }
}
