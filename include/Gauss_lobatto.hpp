#ifndef CARTDG_GAUSS_LOBATTO_HPP_
#define CARTDG_GAUSS_LOBATTO_HPP_

#include "Basis.hpp"

namespace cartdg
{

// Basis based on Gauss-Lobatto quadrature, which is exact for polynomials of degree `2*row_size - 3`
class Gauss_lobatto : public Basis
{
  public:
  Gauss_lobatto (int row_size_arg);
  virtual double node(int i);
  virtual Eigen::VectorXd node_weights();
  virtual Eigen::MatrixXd diff_mat();
  virtual Eigen::MatrixXd boundary();
  virtual Eigen::VectorXd orthogonal(int degree);
  virtual double max_cfl_convective();
  virtual double max_cfl_diffusive();
};

}
#endif
