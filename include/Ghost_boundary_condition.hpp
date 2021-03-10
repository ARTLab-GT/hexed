#ifndef CARTDG_GHOST_BOUNDARY_CONDITION_HPP_
#define CARTDG_GHOST_BOUNDARY_CONDITION_HPP_

#include <Eigen/Dense>
#include <vector>

namespace cartdg
{

class Grid;

class Ghost_boundary_condition
{
  public:
  int i_dim;
  int n_var;
  int n_qpoint;
  bool is_positive_face;
  Eigen::ArrayXXd state; // Stores both states in an order compatible with Flux_kernel
  Eigen::Block<Eigen::ArrayXXd> domain_state(); // indices: (i_qpoint, i_var)
  Eigen::Block<Eigen::ArrayXXd> ghost_state();
  std::vector<double> default_jacobian;
  std::vector<double*> jacobians;
  std::vector<int> elems;

  Ghost_boundary_condition(const Grid& grid, int i_dim_arg, bool is_positive_face_arg);
  virtual ~Ghost_boundary_condition();

  void add_element(int i_elem);
  void add_element(int i_elem, double* jacobian);
  virtual void calc_ghost_state() = 0;
  virtual void print();
};

}
#endif
