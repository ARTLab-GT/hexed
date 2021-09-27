#ifndef CARTDG_GAUSS_LEGENDRE_HPP_
#define CARTDG_GAUSS_LEGENDRE_HPP_

#include "Basis.hpp"

namespace cartdg
{

class Gauss_legendre : public Basis
{
  public:
  Gauss_legendre (int row_size_arg);
  virtual double node(int i);
  virtual Eigen::VectorXd node_weights();
  virtual Eigen::MatrixXd diff_mat();
  virtual Eigen::VectorXd orthogonal(int degree);
};

}
#endif
