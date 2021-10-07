#include <Basis.hpp>

namespace cartdg
{

Basis::Basis(int row_size_arg) : row_size(row_size_arg) {}
Basis::~Basis() {}

void Basis::interpolate(Eigen::VectorXd& values, Eigen::VectorXd& sample, Eigen::VectorXd& write)
{
}

}
