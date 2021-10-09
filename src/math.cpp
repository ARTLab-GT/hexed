#include <math.hpp>

namespace cartdg
{
namespace custom_math
{

void inplace_hcmv(const Eigen::MatrixXd& mat, Eigen::MatrixXd& vec, int i_dim)
{
}

Eigen::VectorXd hypercube_matvec(const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec)
{
  int in_size = 1;
  for (int i_dim = 0, out_size = vec.size(); out_size > 1; ++i_dim)
  {
    #ifdef DEBUG
    if (out_size%mat.cols() != 0)
    {
      const int n {100};
      char buffer [n];
      auto format = "Incompatible matrix shapes in hypercube_matvec (dim %i): cannot divide %i elements into %i rows.";
      snprintf(buffer, n, format, i_dim, out_size, mat.cols());
      throw std::runtime_error(buffer);
    }
    #endif
    in_size *= mat.rows();
    out_size /= mat.cols();
  }
  return vec;
}

}
}
