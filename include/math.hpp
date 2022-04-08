#ifndef CARTDG_MATH_HPP_
#define CARTDG_MATH_HPP_

#include <cmath>
#include <Eigen/Dense>

namespace cartdg::custom_math
{

/*
 * Raises an arbitrary arithmetic type to an integer power. Can return
 * `constexpr`, which `std::pow` is not allowed to do according to the standard
 * (although the GCC implementation can anyway).
 */
template<typename number_t>
constexpr number_t pow(number_t base, int exponent)
{
  number_t result = 1;
  for (int i = 0; i < exponent; ++i) result *= base;
  for (int i = 0; i > exponent; --i) result /= base;
  return result;
}

/*
 * Integer logarithm. If base < 2, returns -1 to indicate failure.
 * Otherwise, if arg < 1, returns 0. In the usual case where neither
 * of the above are true, returns $\ciel{\log_{base}(arg)}$.
 */
constexpr int log(int base, int arg)
{
  if (base <= 1) return -1;
  int result = 0;
  for (int compare = 1; compare < arg; compare *= base) ++result;
  return result;
}

template <typename Func_type>
double root(Func_type func, double init_guess, double atol=1e-10, double init_diff=1e-3)
{
  double guess_prev = init_guess - init_diff;
  double func_prev = func(guess_prev);
  double guess = init_guess;
  do
  {
    double func_curr = func(guess);
    double slope = (func_curr - func_prev)/(guess - guess_prev);
    guess_prev = guess;
    func_prev = func_curr;
    guess -= func_curr/slope;
  }
  while (std::abs(guess - guess_prev) > atol);
  return guess;
}

Eigen::VectorXd hypercube_matvec(const Eigen::MatrixXd&, const Eigen::VectorXd&);
Eigen::VectorXd dimension_matvec(const Eigen::MatrixXd&, const Eigen::VectorXd&, int i_dim);

/*
 * Orthonormalize a basis. Assumes `basis` is invertible. Returns a matrix with the
 * following properties:
 * - Unitary.
 * - Span of columns excluding the `i_dim`th is the same as for `basis`.
 * - Minimizes RMS difference between columns exclucing `i_dim`th of return matrix and `basis`
 *   (sensitive to order).
 * - Inner product of `i_dim`th columns of return matrix and `basis` is positive.
 */
template <int n_dim>
Eigen::Matrix<double, n_dim, n_dim> orthonormal (Eigen::Matrix<double, n_dim, n_dim> basis, int i_dim)
{
  static_assert (n_dim <= 3, "Not implemented for n_dim > 3.");
  if constexpr (n_dim == 1)
  {
    return basis/std::abs(basis(0, 0));
  }
  else
  {
    auto col_i = basis.col(i_dim);
    std::array<int, n_dim - 1> j_col;
    for (int offset = 1; offset < n_dim; ++offset)
    {
      j_col[offset - 1] = (offset + i_dim)%n_dim;
    }
    auto cols {basis(Eigen::all, j_col)}; // all the cols except for `i_dim`th
    cols.array().rowwise() /= cols.array().colwise().norm(); // normalize
    if constexpr (n_dim == 3) // orthonormalize `cols`
    {
      Eigen::Matrix2d sum_diff {{1, -1}, {1, 1}};
      cols = cols*sum_diff;
      cols.array().rowwise() /= cols.array().colwise().norm();
      cols = cols*(sum_diff.transpose()/std::sqrt(2.));
    }
    for (int jc : j_col) // orthogonalize `i_dim`th col wrt `cols` (Gram-Schmidt style)
    {
      col_i -= col_i.dot(basis.col(jc))*basis.col(jc);
    }
    col_i /= col_i.norm();
    return basis;
  }
}

Eigen::MatrixXd orthonormal (Eigen::MatrixXd basis, int i_dim);

}
#endif
