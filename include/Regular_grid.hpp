#ifndef CARTDG_REGULAR_GRID_HPP_
#define CARTDG_REGULAR_GRID_HPP_

#include "Grid.hpp"

namespace cartdg
{

class Regular_grid : public Grid
{
  public:
  Regular_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);
};

}
#endif
