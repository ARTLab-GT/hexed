#ifndef GAUSS_LOBATTO_HPP_
#define GAUSS_LOBATTO_HPP_

#include <Basis.hpp>

class Gauss_lobatto : public Basis
{
  public:
  Gauss_lobatto (int rank_arg);
  virtual double node(int i);
  virtual Eigen::VectorXd node_weights();
  virtual Eigen::MatrixXd diff_mat();
};

#endif