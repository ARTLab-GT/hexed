#include "Equidistant.hpp"

Equidistant::Equidistant(int rank_arg) : Basis(rank_arg) {}

double Equidistant::node(int i)
{
  return 0;
}

double Equidistant::diff_mat(int i_result, int i_operand)
{
  return 0;
}
