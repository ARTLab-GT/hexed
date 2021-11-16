#include <Basis.hpp>

namespace cartdg
{

Basis::Basis(int row_size_arg) : row_size(row_size_arg) {}
Basis::~Basis() {}

Eigen::MatrixXd Basis::interpolate(const Eigen::VectorXd& sample)
{
  Eigen::MatrixXd interp {Eigen::MatrixXd::Ones(sample.size(), row_size)};
  for (int i_node = 0; i_node < row_size; ++i_node)
  {
    for (int j_node = 0; j_node < row_size; ++j_node)
    {
      if (i_node != j_node)
      {
        interp.col(i_node).array() *= (sample - Eigen::VectorXd::Constant(sample.size(), node(j_node))).array()/(node(i_node) - node(j_node));
      }
    }
  }
  return interp;
}

Eigen::MatrixXd Basis::prolong (int i_half)
{
  throw std::runtime_error("Not implemented for base class.");
  return Eigen::MatrixXd::Zero(row_size, row_size);
}

Eigen::MatrixXd Basis::restrict (int i_half)
{
  throw std::runtime_error("Not implemented for base class.");
  return Eigen::MatrixXd::Zero(row_size, row_size);
}

}
