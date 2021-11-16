#ifndef CARTDG_BASIS_HPP_
#define CARTDG_BASIS_HPP_

#include <Eigen/Dense>

namespace cartdg
{

/*
 * Represents a basis of Lagrange (interpolating) polynomials. Interpolates between
 * a set of quadrature points, a.k.a. nodes. The `i`th basis vector is the polynomial with a value
 * of 1 at the `i`th node and 0 at all other nodes.
 */
class Basis
{
  public:
  int row_size;

  Basis(int row_size_arg);
  virtual ~Basis();
  // Return the `i`th interpolation node (i.e. quadrature point)
  // Should be in interval [0, 1].
  virtual double node(int i) = 0;
  // `node_weights(i)` is the quadrature weight associated with `node(i)`
  // Integrates in domain [0, 1].
  virtual Eigen::VectorXd node_weights() = 0;
  // Differentiation matrix. `diff_mat(i, j)` is the derivative of interpolating polynomial `j` evaluated at `node(i)`.
  virtual Eigen::MatrixXd diff_mat() = 0;
  // `boundary()(i, j)` is the `j`th basis polynomial evaluated at `i`. `i` must be in {0, 1}.
  virtual Eigen::MatrixXd boundary() = 0;
  // `orthogonal(degree)(i)` is the Legendre (orthogonal) polynomial evaluated at `node(i)`. Polynomial is
  // scaled so that its norm is 1 with respect to the quadrature rule.
  virtual Eigen::VectorXd orthogonal(int degree) = 0;
  Eigen::MatrixXd interpolate(const Eigen::VectorXd& sample);
  virtual Eigen::MatrixXd prolong  (int i_half);
  virtual Eigen::MatrixXd restrict (int i_half);
};

}
#endif
