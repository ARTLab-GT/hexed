#ifndef HEXED_DERIVATIVE_HPP_
#define HEXED_DERIVATIVE_HPP_

#include <Eigen/Dense>
#include "Basis.hpp"

namespace hexed
{

/*
 * Represents the operator on L_2([0, 1]) which returns the derivative projected
 * onto the space of polynomials of row size `row_size` by discontinuous Galerkin projection.
 */
template <int row_size>
class Derivative
{
  private:
  const Eigen::Matrix<double, row_size, row_size> diff_mat;
  const Eigen::Matrix<double, 2, row_size> boundary;
  const Eigen::Matrix<double, 2, 2> sign {{-1, 0}, {0, 1}};
  const Eigen::Matrix<double, row_size, 1> inv_weights;
  const Eigen::Matrix<double, row_size, 2> lift;

  public:
  Derivative(const Basis& basis)
  : diff_mat {basis.diff_mat()}, boundary {basis.boundary()},
    inv_weights {Eigen::Array<double, row_size, 1>::Constant(1.)/basis.node_weights().array()},
    lift {inv_weights.asDiagonal()*basis.boundary().transpose()*sign}
  {}

  /*
   * `qpoint_vals`: values of a function at the quadrature points
   * `boundary_vals`: values of the function at the boundaries of the interval ({0, 1})
   * returns: values of the derivative at the quadrature points. This derivative is exact
   * if the function is a polynomial of row size `row_size`. It is conservative in the sense
   * that the integral of the return value is equal to the difference between the specified boundary values.
   */
  template<int n_var>
  Eigen::Matrix<double, row_size, n_var> interior_term(const Eigen::Matrix<double, row_size, n_var>& qpoint_vals,
                                                       Eigen::Matrix<double, 2, n_var>& boundary_vals)
  {
    boundary_vals = boundary*qpoint_vals;
    return diff_mat*qpoint_vals - lift*boundary_vals;
  }
  template<int n_var>
  Eigen::Matrix<double, row_size, n_var> boundary_term(const Eigen::Matrix<double, 2, n_var>& boundary_vals)
  {
    return lift*boundary_vals;
  }
  template<int n_var>
  Eigen::Matrix<double, row_size, n_var> operator()(const Eigen::Matrix<double, row_size, n_var>& qpoint_vals,
                                                    const Eigen::Matrix<double, 2, n_var>& boundary_vals)
  {
    return diff_mat*qpoint_vals + lift*(boundary_vals - boundary*qpoint_vals);
  }
  template<int n_var>
  Eigen::Matrix<double, row_size, n_var> operator()(const Eigen::Matrix<double, row_size, n_var>& qpoint_vals)
  {
    return diff_mat*qpoint_vals;
  }
};

}
#endif
