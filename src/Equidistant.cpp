#include <Equidistant.hpp>
#include <assert.hpp>

namespace hexed
{

double Equidistant::min_eig_convection() const
{
  HEXED_ASSERT(false, "Not implemented");
  return 0;
}

double Equidistant::quadratic_safety() const
{
  HEXED_ASSERT(false, "Not implemented");
  return 0;
}

Equidistant::Equidistant(int row_size_arg) : Basis(row_size_arg) {}

double Equidistant::node(int i) const
{
  return (double)i/(row_size - 1);
}

Eigen::VectorXd Equidistant::node_weights() const
{
  throw std::runtime_error("Not implemented.");
  Eigen::VectorXd unused (0);
  return unused;
}

Eigen::MatrixXd Equidistant::diff_mat() const
{
  Eigen::MatrixXd dm (row_size, row_size);
  for (int i_operand = 0; i_operand < row_size; ++i_operand)
  {
    for (int i_result = 0; i_result < row_size; ++i_result)
    {
      double& deriv = dm(i_result, i_operand);
      if (i_result == i_operand)
      {
        deriv = 0.;
        for (int i_node = 0; i_node < row_size; ++i_node)
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
        for (int i_node = 0; i_node < row_size; ++i_node)
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

Eigen::MatrixXd Equidistant::boundary() const
{
  Eigen::MatrixXd b {Eigen::MatrixXd::Zero(2, row_size)};
  b(0, 0) = 1.;
  b(1, row_size - 1) = 1.;
  return b;
}

Eigen::VectorXd Equidistant::orthogonal(int degree) const
{
  HEXED_ASSERT(false, "Not implemented");
  return {};
}

Eigen::MatrixXd Equidistant::filter() const
{
  HEXED_ASSERT(false, "Not implemented");
  return {};
}

Eigen::MatrixXd Equidistant::prolong (int i_half) const
{
  HEXED_ASSERT(false, "Not implemented");
  return {};
}

Eigen::MatrixXd Equidistant::restrict (int i_half) const
{
  HEXED_ASSERT(false, "Not implemented");
  return {};
}

double Equidistant::min_eig_diffusion() const
{
  HEXED_ASSERT(false, "Not implemented");
  return 0;
}

}
