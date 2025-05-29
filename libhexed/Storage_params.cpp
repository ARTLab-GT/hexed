#include <math.hpp>
#include <Storage_params.hpp>

namespace hexed {

int Storage_params::n_qpoint() const {
  return math::pow(row_size, n_dim);
}

int Storage_params::n_face_qpoint() const {
  return math::pow(row_size, n_dim - 1);
}

int Storage_params::n_dof() const {
  return n_qpoint()*n_var;
}

int Storage_params::size() const {
  return n_dof()*n_stage;
}

int Storage_params::n_vertices() const {
  return math::pow(2, n_dim);
}

int Storage_params::n_var_numeric() const {
  return n_var + 3 + n_forcing + n_advection(row_size) + std::max((n_stage - 1)*n_var, n_advection(row_size));
}

int Storage_params::n_dof_numeric() const {
  return n_var_numeric()*n_qpoint();
}

std::vector<int> Storage_params::physical_shape() const {
  std::vector<int> shape(n_dim + 1, row_size);
  shape[0] = n_var;
  return shape;
}

std::vector<int> Storage_params::numerical_shape() const {
  std::vector<int> shape(n_dim + 1, row_size);
  shape[0] = n_var_numeric();
  return shape;
}

bool operator==(Storage_params par0, Storage_params par1) {
  return par0.n_stage == par1.n_stage && par0.n_var == par1.n_var && par0.n_dim == par1.n_dim
         && par0.row_size == par1.row_size && par0.n_forcing == par1.n_forcing;
}

bool operator!=(Storage_params par0, Storage_params par1) {
  return !(par0 == par1);
}

std::string to_string(Storage_params par) {
  return format_str("Storage_params{%istage x %ivar x %idim x %i rows (%iforcing)}",
                    par.n_stage, par.n_var, par.n_dim, par.row_size, par.n_forcing);
}

}
