#ifndef CARTDG_SOLUTION_HPP_
#define CARTDG_SOLUTION_HPP_

#include "Grid.hpp"
#include "Gauss_lobatto.hpp"
#include "kernels/kernel_types.hpp"

namespace cartdg
{

class Initializer;

class Solution
{
  public:
  int n_var;
  int n_dim;
  double base_mesh_size;
  Gauss_lobatto basis;
  std::vector<Grid> grids;

  Solution(int n_var_arg, int n_dim_arg, int rank_arg, double bms);
  virtual ~Solution();

  // functions that access information
  Grid& get_grid(int order_added);
  void visualize(std::string file_prefix);

  // functions that modify the state data
  double update(double cfl_by_stable_cfl=0.7);
  void initialize(Initializer& init);

  // functions that modify the Grid(s)
  void add_block_grid(int ref_level, std::vector<int> lower_corner,
                                     std::vector<int> upper_corner);
  void add_block_grid(int ref_level);
  void add_empty_grid(int ref_level);
  void auto_connect();
  void clear_neighbors();

  protected:
  double refined_mesh_size(int ref_level);

  virtual Local_kernel get_local_kernel();
  virtual Neighbor_kernel get_neighbor_kernel();
  virtual Max_char_speed_kernel get_max_char_speed_kernel();
};

}
#endif
