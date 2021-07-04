#ifndef MATH_HPP_
#define MATH_HPP_

namespace cartdg
{

double root(double (*func)(double arg), double init_guess, double atol=1e-10);

}

#endif
