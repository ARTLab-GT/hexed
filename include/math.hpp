#ifndef MATH_HPP_
#define MATH_HPP_

namespace cartdg
{

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


}

#endif
