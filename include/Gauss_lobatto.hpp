#ifndef HEXED_GAUSS_LOBATTO_HPP_
#define HEXED_GAUSS_LOBATTO_HPP_

#include "Basis.hpp"

namespace hexed
{

// Basis based on Gauss-Lobatto quadrature, which is exact for polynomials of degree `2*row_size - 3`
class Gauss_lobatto : public Basis
{
  public:
  Gauss_lobatto (int row_size_arg);
  virtual double node(int i) const;
  virtual Eigen::VectorXd node_weights() const;
  virtual Eigen::MatrixXd diff_mat() const;
  virtual Eigen::MatrixXd boundary() const;
  virtual Eigen::VectorXd orthogonal(int degree) const;
  virtual double max_cfl_convective() const;
  virtual double max_cfl_diffusive() const;
  virtual double cancellation_convective() const;
  virtual double cancellation_diffusive() const;
};

}
#endif
