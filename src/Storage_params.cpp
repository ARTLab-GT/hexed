#include <math.hpp>
#include <Storage_params.hpp>

namespace cartdg
{

Storage_params::Storage_params(int n_stage, int n_var, int n_dim, int row_size)
: n_stage{n_stage}, n_var{n_var}, row_size{row_size}, n_dim{n_dim},
  n_qpoint{custom_math::pow(row_size, n_dim)}, n_dof{n_qpoint*n_var}, size{n_dof*n_stage}
{}

}
