#include <catch.hpp>
#include <cmath>
#include <math.hpp>

double lin_func(double x)
{
  return 2.*x - 1.;
}
double quad_func(double x)
{
  return (x + .6)*(4.*x - 8.);
}
double exp_func(double x)
{
  return std::exp(x) - 2.;
}

TEST_CASE("root finder")
{
  REQUIRE(cartdg::root(lin_func, -1e7) == Approx(0.5));
  REQUIRE(cartdg::root(quad_func, -0.8) == Approx(-.6));
  REQUIRE(cartdg::root(quad_func, 2.3) == Approx(2.));
  REQUIRE(cartdg::root(exp_func, 0.) == Approx(std::log(2.)));
}
