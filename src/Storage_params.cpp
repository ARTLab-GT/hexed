#include <math.hpp>
#include <Storage_params.hpp>

namespace cartdg
{

unsigned Storage_params::n_qpoint()
{
  return custom_math::pow(row_size, n_dim);
}

unsigned Storage_params::n_dof()
{
  return n_qpoint()*n_var;
}

unsigned Storage_params::size()
{
  return n_dof()*n_stage;
}

}
