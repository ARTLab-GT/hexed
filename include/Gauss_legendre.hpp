#ifndef HEXED_GAUSS_LEGENDRE_HPP_
#define HEXED_GAUSS_LEGENDRE_HPP_

#include "Basis.hpp"

namespace hexed
{

/*! \brief Basis based on [Gauss-Legendre](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature) quadrature.
 * \details The `node`s and `node_weights` of this basis
 * are those of the Gauss-Legendre quadrature,
 * which is exact for polynomials of degree `2*row_size - 1`.
 */
class Gauss_legendre : public Basis
{
  public:
  //! see `hexed::Basis::Basis`.
  Gauss_legendre (int row_size_arg);
  double node(int i) const override;
  Eigen::VectorXd node_weights() const override;
  Eigen::MatrixXd diff_mat() const override;
  Eigen::MatrixXd boundary() const override;
  Eigen::VectorXd orthogonal(int degree) const override;
  Eigen::MatrixXd filter() const override;
  Eigen::MatrixXd prolong  (int i_half) const override;
  Eigen::MatrixXd restrict (int i_half) const override;
  double max_cfl_convective() const override;
  double max_cfl_diffusive() const override;
  double cancellation_convective() const override;
  double cancellation_diffusive() const override;
};

}
#endif
