#ifndef EQUIDISTANT_HPP_
#define EQUIDISTANT_HPP_

#include "Basis.hpp"

class Equidistant : public Basis
{
  public:
  Equidistant(int rank_arg);
  virtual double node(int i);
  virtual double diff_mat(int i_result, int i_operand);
};

#endif
