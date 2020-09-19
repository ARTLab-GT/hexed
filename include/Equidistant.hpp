#ifndef EQUIDISTANT_HPP_
#define EQUIDISTANT_HPP_

#include "Basis.hpp"

class Equidistant : public Basis
{
  public:
  Equidistant(int rank_arg);
  virtual double node(int i);
  Eigen::MatrixXd diff_mat();
};

#endif
