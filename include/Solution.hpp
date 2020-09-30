#ifndef SOLUTION_HPP_
#define SOLUTION_HPP_

#include <Grid.hpp>
#include <Gauss_lobatto.hpp>
#include <kernels/kernel_types.hpp>

class Initializer;

class Solution
{
  public:
  int n_var;
  int n_dim;
  double base_mesh_size;

  Solution(int n_var_arg, int n_dim_arg, int rank_arg, double bms);
  virtual ~Solution();

  Grid& get_grid(int order_added);
  double update(double cfl_by_stable_cfl=0.7);

  void initialize(Initializer& init);
  void add_block_grid(int ref_level, std::vector<int> lower_corner,
                                     std::vector<int> upper_corner);
  void add_block_grid(int ref_level);
  void auto_connect();

  protected:
  Gauss_lobatto basis;
  std::vector<Grid> grids;
  virtual Local_kernel get_local_kernel();
  virtual Neighbor_kernel get_neighbor_kernel();
  virtual Max_char_speed_kernel get_max_char_speed_kernel();
};

#endif
