#ifndef CARTDG_SOLUTION_HPP_
#define CARTDG_SOLUTION_HPP_

#include "Regular_grid.hpp"
#include "Deformed_grid.hpp"
#include "Gauss_legendre.hpp"
#include "Kernel_settings.hpp"
#include "Spacetime_func.hpp"

namespace cartdg
{

class Solution
{
  public:
  int n_var;
  int n_dim;
  double base_mesh_size;
  Gauss_legendre basis;
  Kernel_settings kernel_settings;
  std::vector<Regular_grid> reg_grids;
  std::vector<Deformed_grid> def_grids;

  Solution(int n_var_arg, int n_dim_arg, int row_size_arg, double bms);
  virtual ~Solution();

  // functions that access information
  void visualize_field(std::string name);
  void visualize_surface(std::string name);
  std::vector<double> integral();
  std::vector<double> integral(Domain_func& integrand);
  std::vector<double> surface_integral(Surface_func& integrand);
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

  protected:
  double time;
  double refined_mesh_size(int ref_level);
};

}
#endif
