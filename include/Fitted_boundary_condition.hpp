#ifndef CARTDG_FITTED_BOUNDARY_CONDITION_HPP_
#define CARTDG_FITTED_BOUNDARY_CONDITION_HPP_

#include <Eigen/Dense>

namespace cartdg
{

class Grid;

class Fitted_boundary_condition
{
  public:
  int i_dim;
  int n_var;
  int n_qpoint;
  bool is_positive_face;
  Eigen::ArrayXXd state; // Stores both states in an order compatible with Flux_kernel
  Eigen::Block<Eigen::ArrayXXd> domain_state(); // indices: (i_qpoint, i_var)
  Eigen::Block<Eigen::ArrayXXd> ghost_state();
  std::vector<int> elems;

  Fitted_boundary_condition(int n_var_arg, int n_qpoint_arg, int i_dim_arg,
                            bool is_positive_face_arg);
  Fitted_boundary_condition(Grid& grid, int i_dim_arg, bool is_positive_face_arg);
  virtual ~Fitted_boundary_condition();

  virtual void calc_ghost_state() = 0;
};

}
#endif
