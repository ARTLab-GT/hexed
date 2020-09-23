#ifndef SOLUTION_HPP_
#define SOLUTION_HPP_

#include "Grid.hpp"
#include "Equidistant.hpp"

class Solution
{
  public:
  int n_var;
  int n_dim;
  double base_mesh_size;

  Solution(int n_var_arg, int n_dim_arg, int rank_arg, double bms);
  virtual ~Solution();

  void add_block_grid(int ref_level, std::vector<int> lower_corner,
                                     std::vector<int> upper_corner);
  Grid& get_grid(int order_added);
  void update();

  protected:
  Equidistant basis;
  std::vector<Grid> grids;
};

#endif
