#ifndef HEXED_GAUSS_LEGENDRE_HPP_
#define HEXED_GAUSS_LEGENDRE_HPP_

#include "Basis.hpp"

namespace hexed
{

// Basis based on Gauss-Legendre quadrature, which is exact for polynomials of degree `2*row_size - 1`
class Gauss_legendre : public Basis
{
  public:
  Gauss_legendre (int row_size_arg);
  virtual double node(int i) const;
  virtual Eigen::VectorXd node_weights() const;
  virtual Eigen::MatrixXd diff_mat() const;
  virtual Eigen::MatrixXd boundary() const;
  virtual Eigen::VectorXd orthogonal(int degree) const;
  virtual Eigen::MatrixXd prolong  (int i_half) const;
  virtual Eigen::MatrixXd restrict (int i_half) const;
  virtual double max_cfl_convective() const;
  virtual double max_cfl_diffusive() const;
};

}
#endif
