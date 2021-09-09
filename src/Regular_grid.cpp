#include <Regular_grid.hpp>

namespace cartdg
{

Regular_grid::Regular_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: Grid (n_var_arg, n_dim_arg, n_elem_arg, mesh_size_arg, basis_arg)
{}

}
