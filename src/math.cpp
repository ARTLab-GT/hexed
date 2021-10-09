#include <math.hpp>

namespace cartdg
{
namespace custom_math
{

Eigen::VectorXd hcmv(const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec, int prod_size)
{
  const int vec_row_size = vec.size()/mat.cols();
  if (vec_row_size == 1) return mat*vec;
  else
  {
    const int row_size = prod_size/mat.rows();
    Eigen::VectorXd fact {mat.cols()*row_size};
    for (int i_row = 0; i_row < mat.cols(); ++i_row)
    {
      fact(Eigen::seqN(i_row*row_size, row_size)) = hcmv(mat, vec(Eigen::seqN(i_row*vec_row_size, vec_row_size)), row_size);
    }
    Eigen::VectorXd prod {prod_size};
    for (int i_col = 0; i_col < row_size; ++i_col)
    {
      Eigen::VectorXd col = fact(Eigen::seqN(i_col, mat.cols(), row_size));
      prod(Eigen::seqN(i_col, mat.rows(), row_size)) = mat*col;
    }
    return prod;
  }
}

Eigen::VectorXd hypercube_matvec(const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec)
{
  int out_size = 1;
  int i_dim = 0;
  for (int in_size = vec.size(); in_size > 1; ++i_dim)
  {
    #ifdef DEBUG
    if (in_size%mat.cols() != 0)
    {
      const int n {100};
      char buffer [n];
      auto format = "Incompatible matrix shapes in hypercube_matvec (dim %i): cannot divide %i elements into %i rows.";
      snprintf(buffer, n, format, i_dim, in_size, mat.cols());
      throw std::runtime_error(buffer);
    }
    #endif
    out_size *= mat.rows();
    in_size /= mat.cols();
  }
  return hcmv(mat, vec, out_size);
}

}
}
