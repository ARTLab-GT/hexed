#ifndef STORAGE_PARAMS_HPP_
#define STORAGE_PARAMS_HPP_

namespace cartdg
{

struct Storage_params
{
  const int n_stage;
  const int n_var;
  const int row_size;
  const int n_dim;
  const int n_qpoint;
  const int n_dof;
  const int size;

  Storage_params(int n_stage, int n_var, int n_dim, int row_size);
};

}

#endif
