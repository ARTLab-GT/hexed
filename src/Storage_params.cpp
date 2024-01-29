#include <math.hpp>
#include <Storage_params.hpp>

namespace hexed
{

int Storage_params::n_qpoint() const
{
  return math::pow(row_size, n_dim);
}

int Storage_params::n_face_qpoint() const
{
  return math::pow(row_size, n_dim - 1);
}

int Storage_params::n_dof() const
{
  return n_qpoint()*n_var;
}

int Storage_params::size() const
{
  return n_dof()*n_stage;
}

int Storage_params::n_vertices() const
{
  return math::pow(2, n_dim);
}

int Storage_params::n_var_numeric() const
{
  return n_var + 3 + n_forcing + n_advection(row_size) + std::max((n_stage - 1)*n_var, n_advection(row_size));
}

int Storage_params::n_dof_numeric() const
{
  return n_var_numeric()*n_qpoint();
}

}
