#ifndef HEXED_EQUIDISTANT_HPP_
#define HEXED_EQUIDISTANT_HPP_

#include "Basis.hpp"

namespace hexed
{

class Equidistant : public Basis
{
  public:
  Equidistant(int row_size_arg);
  virtual double node(int i) const;
  virtual Eigen::MatrixXd diff_mat() const;
  virtual Eigen::VectorXd node_weights() const;
  virtual Eigen::MatrixXd boundary() const;
  virtual Eigen::VectorXd orthogonal(int degree) const;
};

}
#endif
