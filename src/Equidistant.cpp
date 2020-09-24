#include <Equidistant.hpp>

Equidistant::Equidistant(int rank_arg) : Basis(rank_arg) {}

double Equidistant::node(int i)
{
  return (double)i/(rank - 1);
}

Eigen::MatrixXd Equidistant::diff_mat()
{
  Eigen::MatrixXd dm (rank, rank);
  for (int i_operand = 0; i_operand < rank; ++i_operand)
  {
    for (int i_result = 0; i_result < rank; ++i_result)
    {
      double& deriv = dm(i_result, i_operand);
      if (i_result == i_operand)
      {
        deriv = 0.;
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
    }
  }
  return dm;
}

Eigen::VectorXd Equidistant::node_weights()
{
  throw "Not implemented.";
  Eigen::VectorXd unused (0);
  return unused;
}
