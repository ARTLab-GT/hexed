#include <math.hpp>

namespace cartdg
{
namespace custom_math
{

Eigen::VectorXd hypercube_matvec(const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec)
{
  #if DEBUG
  if (vec.size()%mat.cols() != 0)
  {
    const int n {100};
    char buffer [n];
    auto format = "Incompatible matrix shapes in hypercube_matvec: cannot divide row of %i elements into %i sub-rows.";
    snprintf(buffer, n, format, vec.size(), mat.cols());
    throw std::runtime_error(buffer);
  }
  #endif

  const int vec_row_size = vec.size()/mat.cols();
  if (vec_row_size == 1) return mat*vec;
  else
  {
    int prod_size = 1;
    for (int vec_size = vec.size(); vec_size > 1;)
    {
      prod_size *= mat.rows();
      vec_size /= mat.cols();
    }
    const int row_size = prod_size/mat.rows();
    Eigen::VectorXd fact {mat.cols()*row_size};
    for (int i_row = 0; i_row < mat.cols(); ++i_row)
    {
      fact(Eigen::seqN(i_row*row_size, row_size)) = hypercube_matvec(mat, vec(Eigen::seqN(i_row*vec_row_size, vec_row_size)));
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

Eigen::VectorXd dimension_matvec(const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec, int i_dim)
{
  #if DEBUG
  int n_rows {pow(int(mat.cols()), i_dim + 1)};
  if (vec.size()%n_rows)
  {
    const int n {100};
    char buffer [n];
    auto format = "Incompatible matrix shapes in dimension_matvec: cannot divide %i elements into %i^%i = %i rows.";
    snprintf(buffer, n, format, vec.size(), mat.cols(), i_dim + 1, n_rows);
    throw std::runtime_error(buffer);
  }
  #endif
  return Eigen::VectorXd {};
}

}
}
