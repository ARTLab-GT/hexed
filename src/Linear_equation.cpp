#include <Linear_equation.hpp>

namespace hexed
{

double Linear_equation::norm(int input) {return std::sqrt(inner(input, input));}

Dense_equation::Dense_equation(Mat<dyn, dyn> matrix, Mat<> rhs, int nv, Mat<> guess)
: _mat{matrix}, _vecs{Mat<dyn, dyn>::Zero(matrix.cols(), nv)}
{
  HEXED_ASSERT(matrix.rows() == matrix.cols() && matrix.rows() == rhs.size() && matrix.rows() == guess.size(),
               "size of input matrices and vectors must match");
  HEXED_ASSERT(nv >= 2, "must at least provide storage for initial guess and residual");
  _vecs(all, 0) = guess;
  _vecs(all, 1) = rhs - matrix*guess;
}

Dense_equation::Dense_equation(Mat<dyn, dyn> matrix, Mat<> rhs, int nv)
: Dense_equation(matrix, rhs, nv, Mat<>::Zero(rhs.size()))
{}

int Dense_equation::n_vecs() {return _vecs.cols();}

void Dense_equation::scale(int output, int input, double scalar)
{
  vec(output) = scalar*vec(input);
}

void Dense_equation::add(int output, double coef0, int vec0, double coef1, int vec1)
{
  vec(output) = coef0*vec(vec0) + coef1*vec(vec1);
}

double Dense_equation::inner(int input0, int input1)
{
  return vec(input0).dot(vec(input1));
}

void Dense_equation::matvec(int output, int input)
{
  vec(output) = _mat*vec(input);
}

}
