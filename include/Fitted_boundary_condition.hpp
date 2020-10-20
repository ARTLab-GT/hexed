#ifndef CARTDG_FITTED_BOUNDARY_CONDITION_HPP_
#define CARTDG_FITTED_BOUNDARY_CONDITION_HPP_

#include <Eigen/Dense>

namespace cartdg
{

class Fitted_boundary_condition
{
  public:
  Eigen::ArrayXXd state;
  Eigen::ArrayXXd ghost_state;
  int i_dim;
  bool is_positive_face;
  std::vector<int> elems;

  Fitted_boundary_condition(int n_var, int n_qpoint, int i_dim_arg, bool is_positive_face_arg);
  virtual ~Fitted_boundary_condition();

  virtual void calc_ghost_state() = 0;
};

}
#endif
