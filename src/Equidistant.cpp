#include "Equidistant.hpp"

Equidistant::Equidistant(int rank_arg) : Basis(rank_arg) {}

double Equidistant::node(int i)
{
  return (double)i/(rank - 1);
}

double Equidistant::diff_mat(int i_result, int i_operand)
{
  double deriv;
  if (i_result == i_operand)
  {
    deriv = 0;
    for (int i_node = 0; i_node < rank; ++i_node)
    {
      if (i_node != i_operand)
      {
        deriv += 1./(node(i_result) - node(i_node));
      }
    }
  }
  else
  {
    deriv = 1;
    for (int i_node = 0; i_node < rank; ++i_node)
    {
      if (i_node != i_operand)
      {
        deriv /= node(i_operand) - node(i_node);
        if (i_node != i_result)
        {
          deriv *= node(i_result) - node(i_node);
        }
      }
    }
  }
  return deriv;
}
