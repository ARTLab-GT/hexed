#ifndef MATH_HPP_
#define MATH_HPP_

#include <cmath>
#include <Eigen/Dense>

namespace cartdg
{
namespace custom_math
{

template<typename number_t>
constexpr number_t pow(number_t base, int exponent)
{
  number_t result = 1;
  for (int i = 0; i < exponent; ++i) result *= base;
  for (int i = 0; i > exponent; --i) result /= base;
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
Eigen::Matrix<double, n_dim, n_dim> orthonormal (const Eigen::Matrix<double, n_dim, n_dim>& basis, int i_dim)
{
  static_assert (n_dim <= 3, "Not implemented for n_dim > 3.");
  if constexpr (n_dim == 1)
  {
    return basis/std::abs(basis(0, 0));
  }
  else
  {
    Eigen::Matrix<double, n_dim, n_dim> orth {basis};
    auto col_i = orth.col(i_dim);
    for (int offset = 1; offset < n_dim; ++offset)
    {
      auto col_j = orth.col((offset + i_dim)%n_dim);
      col_j /= col_j.norm();
    }
    for (int offset = 1; offset < n_dim; ++offset)
    {
      auto col_j = orth.col((offset + i_dim)%n_dim);
      col_i -= col_i.dot(col_j)*col_j;
    }
    col_i /= col_i.norm();
    return orth;
  }
}

}
}
#endif
