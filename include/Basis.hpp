#ifndef BASIS_HPP_
#define BASIS_HPP_

class Basis
{
  public:
  int rank;

  Basis(int rank_arg);
  virtual ~Basis();
  virtual double node(int i) = 0;
  virtual double diff_mat(int i_result, int i_operand) = 0;
};

#endif
