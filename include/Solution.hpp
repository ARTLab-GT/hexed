#ifndef CARTDG_SOLUTION_HPP_
#define CARTDG_SOLUTION_HPP_

#include "Grid.hpp"
#include "Deformed_grid.hpp"
#include "Gauss_lobatto.hpp"
#include "kernels/kernel_types.hpp"
#include "kernels/Kernel_settings.hpp"
#include "Spacetime_func.hpp"

namespace cartdg
{

class Solution
{
  public:
  int n_var;
  int n_dim;
  double base_mesh_size;
  Gauss_lobatto basis;
  Kernel_settings kernel_settings;
  std::vector<Grid> grids;
  std::vector<Deformed_grid> def_grids;

  Solution(int n_var_arg, int n_dim_arg, int rank_arg, double bms);
  virtual ~Solution();

  // functions that access information
  Grid& get_grid(int order_added);
  void visualize(std::string file_prefix);
  std::vector<double> integral();
  std::vector<double> integral(Domain_func& integrand);
  std::vector<Grid*> all_grids();

  // functions that modify the state data
  double update(double cfl_by_stable_cfl=0.7);
  void initialize(Spacetime_func& init_cond);

  // functions that modify the Grid(s)
  void add_block_grid(int ref_level, std::vector<int> lower_corner,
                                     std::vector<int> upper_corner);
  void add_block_grid(int ref_level);
  void add_empty_grid(int ref_level);
  void add_deformed_grid(int ref_level);
  void auto_connect();
  void clear_neighbors();

  protected:
  double refined_mesh_size(int ref_level);

  virtual Local_kernel get_local_kernel();
  virtual Local_deformed_kernel get_local_deformed_kernel();
  virtual Neighbor_kernel get_neighbor_kernel();
  virtual Max_char_speed_kernel get_max_char_speed_kernel();
  virtual Fbc_kernel get_fbc_kernel();
};

}
#endif
