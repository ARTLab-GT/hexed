#include <catch.hpp>
#include <cmath>
#include <math.hpp>

TEST_CASE("static_pow")
{
  {
    constexpr int result = cartdg::custom_math::pow(2, 0);
    REQUIRE(result == 1);
  }
  {
    constexpr int result = cartdg::custom_math::pow(2, 1);
    REQUIRE(result == 2);
  }
  {
    constexpr int result = cartdg::custom_math::pow(2, -1);
    REQUIRE(result == 0);
  }
  {
    constexpr int result = cartdg::custom_math::pow(3, 4);
    REQUIRE(result == 81);
  }
  {
    constexpr double result = cartdg::custom_math::pow(1.5, 2);
    REQUIRE(result == 2.25);
  }
  {
    constexpr double result = cartdg::custom_math::pow(2., -1);
    REQUIRE(result == 0.5);
  }
}

double lin_func(double x)
{
  return 2.*x - 1.;
}
double quad_func(double x)
{
  return (x + .6)*(4.*x - 8.);
}

TEST_CASE("root finder")
{
  REQUIRE(cartdg::custom_math::root(lin_func, -1e7) == Approx(0.5));
  REQUIRE(cartdg::custom_math::root(quad_func, -0.8) == Approx(-.6));
  REQUIRE(cartdg::custom_math::root(quad_func, 2.3) == Approx(2.));
  REQUIRE(cartdg::custom_math::root([](double x){return std::exp(x) - 2.;}, 0.)
          == Approx(std::log(2.)));
}
