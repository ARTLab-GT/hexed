#ifndef CARTDG_EQUIDISTANT_HPP_
#define CARTDG_EQUIDISTANT_HPP_

#include "Basis.hpp"

namespace cartdg
{

class Equidistant : public Basis
{
  public:
  Equidistant(int rank_arg);
  virtual double node(int i);
  virtual Eigen::MatrixXd diff_mat();
  virtual Eigen::VectorXd node_weights();
  virtual Eigen::VectorXd orthogonal(int degree);
};

}
#endif
