#ifndef STORAGE_PARAMS_HPP_
#define STORAGE_PARAMS_HPP_

namespace cartdg
{

class Storage_params
{
  public:
  int n_stage;
  int n_var;
  int n_dim;
  int row_size;

  int n_qpoint();
  int n_dof();
  int size();
  int n_vertices();
};

}

#endif
