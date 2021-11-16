#ifndef CARTDG_GAUSS_LOBATTO_HPP_
#define CARTDG_GAUSS_LOBATTO_HPP_

#include "Basis.hpp"

namespace cartdg
{

class Gauss_lobatto : public Basis
{
  public:
  Gauss_lobatto (int row_size_arg);
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
