#ifndef CARTDG_DERIVATIVE_HPP_
#define CARTDG_DERIVATIVE_HPP_

#include <Eigen/Dense>
#include "Basis.hpp"

namespace cartdg
{

template <int row_size>
class Derivative
{
  private:
  const Eigen::Matrix<double, row_size, row_size> diff_mat;
  const Eigen::Matrix<double, 2, row_size> boundary;
  const Eigen::Matrix<double, 2, 2> sign {{1, 0}, {0, -1}};
  const Eigen::Matrix<double, row_size, 1> inv_weights;
  const Eigen::Matrix<double, row_size, 2> lift;

  public:
  Derivative(Basis& basis)
  : diff_mat {basis.diff_mat()}, boundary {basis.boundary()},
    inv_weights {Eigen::Array<double, row_size, 1>::Constant(1.)/basis.node_weights().array()},
    lift {inv_weights.asDiagonal()*basis.boundary().transpose()*sign}
  {}

  template<int n_var>
  Eigen::Matrix<double, row_size, n_var> operator()(Eigen::Matrix<double, row_size, n_var>& qpoint_vals,
                                                    Eigen::Matrix<double, 2, n_var>& boundary_vals)
  {
    return Eigen::Matrix<double, row_size, n_var>::Zero();
  }
};

}
#endif
