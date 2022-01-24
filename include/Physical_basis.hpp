#ifndef PHYSICAL_BASIS_HPP_
#define PHYSICAL_BASIS_HPP_

#include <vector>
#include <Eigen/Dense>

namespace cartdg
{

/*
 * Represents a polynomial basis in physical space (as opposed to reference
 * space). The basis functions can be evaluated at a set of quadrature
 * points. It spans the space of polynomials where the total powers of
 * all variables in each term is less than `row_size`. That is, the same
 * space traditionally used for simplex elements. For example,
 * x_0^(`row_size` - 1) and x_0^(`row_size` - 3)*x_1^(2) are in the space
 * but x_0^('row_size` - 1)*x_1^(`row_size` - 1) is not. This basis will
 * be useful for handling degenerate elements.
 */
class Physical_basis
{
  int n_dim;
  int row_size;
  int n_qpoint;
  std::vector<double> pos;
  std::vector<double> bound_box_center; // coordinates of center of bounding box of points
  std::vector<double> bound_box_size; // side length in each dimension
  std::vector<std::vector<double>> indices;

  public:
  /*
   * `qpoint_pos` represents position of quadrature points. Must satisfy `qpoint_pos.size() == n_dim*pow(row_size, n_dim)`.
   * `n_dim` must be <= 3.
   * Any qpoint distribution is allowed, but some might result in a singular basis.
   * If efficiency becomes a problem, could pass `pos` by reference
   */
  Physical_basis(int n_dim, int row_size, std::vector<double> pos);
  int size(); // return the size of the basis, i.e. the number of basis polynomials (not the same as n_qpoint).
  /*
   * Evaluate the `i_basis`th polynomial at the `i_qpoint`th quadrature point.
   * Args must satisfy `0 <= i_qpoint < n_qpoint; 0 <= i_basis < size()`.
   */
  double evaluate(int i_qpoint, int i_basis);
  /*
   * Compute the projection of `polys` onto the span of the physical basis.
   * Each column of `polys` represents a polynomial to be projected, evaluated at
   * the quadrature points. Projection is orthogonal with respect to quadrature
   * using `weights` . Requires `polys.rows()`, `weights.size()`, and `pos.size()`
   * to be equal. Returns projection evaluated at quadrature points.
   */
  Eigen::MatrixXd projection(Eigen::MatrixXd polys, Eigen::VectorXd weights);
};

}
#endif
