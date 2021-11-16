#ifndef CARTDG_BASIS_HPP_
#define CARTDG_BASIS_HPP_

#include <Eigen/Dense>

namespace cartdg
{

class Basis
{
  public:
  int row_size;

  Basis(int row_size_arg);
  virtual ~Basis();
  virtual double node(int i) = 0;
  virtual Eigen::VectorXd node_weights() = 0;
  virtual Eigen::MatrixXd diff_mat() = 0;
  virtual Eigen::MatrixXd boundary() = 0;
  virtual Eigen::VectorXd orthogonal(int degree) = 0;
  Eigen::MatrixXd interpolate(const Eigen::VectorXd& sample);
  virtual Eigen::MatrixXd prolong  (int i_half);
  virtual Eigen::MatrixXd restrict (int i_half);
};

}
#endif
