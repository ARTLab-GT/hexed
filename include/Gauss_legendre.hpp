#ifndef CARTDG_GAUSS_LEGENDRE_HPP_
#define CARTDG_GAUSS_LEGENDRE_HPP_

#include "Basis.hpp"

namespace cartdg
{

// Basis based on Gauss-Legendre quadrature, which is exact for polynomials of degree `2*row_size - 1`
class Gauss_legendre : public Basis
{
  public:
  Gauss_legendre (int row_size_arg);
  virtual double node(int i);
  virtual Eigen::VectorXd node_weights();
  virtual Eigen::MatrixXd diff_mat();
  virtual Eigen::MatrixXd boundary();
  virtual Eigen::VectorXd orthogonal(int degree);
  virtual Eigen::MatrixXd prolong  (int i_half);
  virtual Eigen::MatrixXd restrict (int i_half);
};

}
#endif
