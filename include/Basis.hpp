#ifndef CARTDG_BASIS_HPP_
#define CARTDG_BASIS_HPP_

#include <Eigen/Dense>

namespace cartdg
{

class Basis
{
  public:
  int rank;

  Basis(int rank_arg);
  virtual ~Basis();
  virtual double node(int i) = 0;
  virtual Eigen::VectorXd node_weights() = 0;
  virtual Eigen::MatrixXd diff_mat() = 0;
};

}
#endif
