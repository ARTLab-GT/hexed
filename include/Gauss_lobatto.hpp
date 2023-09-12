#ifndef HEXED_GAUSS_LOBATTO_HPP_
#define HEXED_GAUSS_LOBATTO_HPP_

#include "Basis.hpp"

namespace hexed
{

/*! \brief Basis based on [Gauss-Lobatto](https://mathworld.wolfram.com/LobattoQuadrature.html) quadrature.
 * \details The `node`s and `node_weights` of this basis
 * are those of the Gauss-Legendre quadrature,
 * which is exact for polynomials of degree `2*row_size - 3`.
 */
class Gauss_lobatto : public Basis
{
  protected:
  double min_eig_convection() const override;
  inline double quadratic_safety() const override {return 0.7;}
  public:
  //! see `hexed::Basis::Basis`.
  Gauss_lobatto (int row_size_arg);
  double node(int i) const override;
  Eigen::VectorXd node_weights() const override;
  Eigen::MatrixXd diff_mat() const override;
  Eigen::MatrixXd boundary() const override;
  Eigen::VectorXd orthogonal(int degree) const override;
  Eigen::MatrixXd filter() const override;
  Eigen::MatrixXd prolong(int i_half) const override;
  Eigen::MatrixXd restrict(int i_half) const override;
  double min_eig_diffusion() const override;
};

}
#endif
