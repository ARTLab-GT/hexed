#include <iostream>

#include <catch2/catch.hpp>
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

TEST_CASE("hypercube_matvec")
{
  auto hcmv {cartdg::custom_math::hypercube_matvec};
  #ifdef DEBUG
  SECTION("multiplying incompatible shapes throws")
  {
    {
      Eigen::MatrixXd mat {Eigen::MatrixXd::Identity(6, 5)};
      Eigen::VectorXd vec {Eigen::VectorXd::Ones(4)};
      REQUIRE_THROWS(hcmv(mat, vec));
    }
    {
      Eigen::MatrixXd mat {Eigen::MatrixXd::Identity(2, 3)};
      Eigen::VectorXd vec {Eigen::VectorXd::Ones(28)};
      REQUIRE_THROWS(hcmv(mat, vec));
    }
    {
      Eigen::MatrixXd mat {Eigen::MatrixXd::Identity(2, 3)};
      Eigen::VectorXd vec {Eigen::VectorXd::Ones(36)};
      REQUIRE_THROWS(hcmv(mat, vec));
    }
  }
  #endif
  SECTION("correct values")
  {
    Eigen::MatrixXd mat {{0.5, 0.5, 0.}, {0., 0.5, 0.5}};
    Eigen::VectorXd vec {Eigen::VectorXd::LinSpaced(27, 0, 26)};
    auto prod = hcmv(mat, vec);
    REQUIRE(prod.size() == 8);
    Eigen::VectorXd correct {{6.5, 7.5, 9.5, 10.5, 15.5, 16.5, 18.5, 19.5}};
    REQUIRE(prod == correct);
  }
}
