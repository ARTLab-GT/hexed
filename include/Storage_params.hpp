#ifndef STORAGE_PARAMS_HPP_
#define STORAGE_PARAMS_HPP_

namespace cartdg
{

class Storage_params
{
  public:
  unsigned n_stage;
  unsigned n_var;
  unsigned n_dim;
  unsigned row_size;

  unsigned n_qpoint();
  unsigned n_dof();
  unsigned size();
};

}

#endif
