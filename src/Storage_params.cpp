#include <math.hpp>
#include <Storage_params.hpp>

namespace hexed
{

int Storage_params::n_qpoint()
{
  return custom_math::pow(row_size, n_dim);
}

int Storage_params::n_dof()
{
  return n_qpoint()*n_var;
}

int Storage_params::size()
{
  return n_dof()*n_stage;
}

int Storage_params::n_vertices()
{
  return custom_math::pow(2, n_dim);
}

}
