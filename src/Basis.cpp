#include <Basis.hpp>

namespace cartdg
{

Basis::Basis(int row_size_arg) : row_size(row_size_arg) {}
Basis::~Basis() {}

Eigen::MatrixXd Basis::interpolate(Eigen::VectorXd& sample)
{
  return Eigen::MatrixXd::Identity(sample.size(), row_size);
}

}
