#include <iostream>

#include <catch2/catch.hpp>
#include <cmath>
#include <hexed/math.hpp>

// yo, compile time tests! Here we test constexpr `pow` and `log`
static_assert (hexed::custom_math::pow(2, 0) == 1);
static_assert (hexed::custom_math::pow(2, 1) == 2);
static_assert (hexed::custom_math::pow(2, -1) == 0);
static_assert (hexed::custom_math::pow(3, 4) == 81);
static_assert (hexed::custom_math::pow(1.5, 2) == 2.25);
static_assert (hexed::custom_math::pow(2., -1) == 0.5);
static_assert (hexed::custom_math::log(2, -1) == 0);
static_assert (hexed::custom_math::log(2, 1) == 0);
static_assert (hexed::custom_math::log(2, 2) == 1);
static_assert (hexed::custom_math::log(2, 16) == 4);
static_assert (hexed::custom_math::log(2, 15) == 4);
static_assert (hexed::custom_math::log(3, 27) == 3);
static_assert (hexed::custom_math::log(1, 27) == -1);
static_assert (hexed::custom_math::log(-1, 27) == -1);

double lin_func(double x)
{
  return 2.*x - 1.;
}
double quad_func(double x)
{
  return (x + .6)*(4.*x - 8.);
}

TEST_CASE("broyden root finder")
{
  REQUIRE(hexed::custom_math::broyden(lin_func, -1e7) == Approx(0.5));
  REQUIRE(hexed::custom_math::broyden(quad_func, -0.8) == Approx(-.6));
  REQUIRE(hexed::custom_math::broyden(quad_func, 2.3) == Approx(2.));
  REQUIRE(hexed::custom_math::broyden([](double x){return std::exp(x) - 2.;}, 0.)
          == Approx(std::log(2.)));
}

TEST_CASE("bisection root finder")
{
  REQUIRE(hexed::custom_math::bisection(lin_func , {  0,   2}) == Approx(0.5));
  REQUIRE(hexed::custom_math::bisection(quad_func, { -1,   0}) == Approx(-.6));
  REQUIRE(hexed::custom_math::bisection(quad_func, {0.1, 3.4}) == Approx(2.));
  REQUIRE(hexed::custom_math::bisection([](double x){return std::exp(x) - 2.;}, {0, 1})
          == Approx(std::log(2.)));
}

TEST_CASE("hypercube_matvec")
{
  auto hcmv {hexed::custom_math::hypercube_matvec};
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

TEST_CASE("dimension matvec")
{
  auto dmv {hexed::custom_math::dimension_matvec};
  #ifdef DEBUG
  SECTION("multiplying incompatible shapes throws")
  {
    {
      Eigen::MatrixXd mat {Eigen::MatrixXd::Identity(7, 7)};
      Eigen::VectorXd vec {Eigen::VectorXd::Zero(10)};
      REQUIRE_THROWS(dmv(mat, vec, 0));
    }
    {
      Eigen::MatrixXd mat {Eigen::MatrixXd::Identity(4, 5)};
      Eigen::VectorXd vec {Eigen::VectorXd::Zero(750)};
      dmv(mat, vec, 2);
      REQUIRE_THROWS(dmv(mat, vec, 3));
    }
  }
  #endif
  Eigen::VectorXd vec {Eigen::VectorXd::LinSpaced(8, 0, 7)};
  Eigen::MatrixXd mat {{0, 1}, {1, 0}, {0.5, 0.5}};
  auto prod = dmv(mat, vec, 1);
  REQUIRE(prod.size() == 12);
  Eigen::VectorXd correct {{2, 3, 0, 1, 1, 2, 6, 7, 4, 5, 5, 6}};
  REQUIRE(prod == correct);
}

TEST_CASE("orthonormal")
{
  SECTION("1D")
  {
    Eigen::Matrix<double, 1, 1> basis {-1.2};
    REQUIRE(hexed::custom_math::orthonormal(basis, 0)(0, 0) == Approx(-1.));
    basis(0, 0) = 0.1;
    REQUIRE(hexed::custom_math::orthonormal(basis, 0)(0, 0) == Approx(1.));
  }
  SECTION("2D")
  {
    Eigen::Matrix2d basis;
    basis << -2., 0.3,
              0., 0.4;
    Eigen::Matrix2d backup = basis*1.;
    Eigen::Matrix2d correct;
    correct << -1., 0.,
                0., 1.;
    REQUIRE((hexed::custom_math::orthonormal(basis, 1) - correct).norm() == Approx(0.).scale(1.));
    REQUIRE((backup - basis).norm() == Approx(0.).scale(1.)); // make sure we don't modify `basis`
    correct << -0.8, 0.6,
                0.6, 0.8;
    REQUIRE((hexed::custom_math::orthonormal(basis, 0) - correct).norm() == Approx(0.).scale(1.));
  }
  SECTION("3D")
  {
    Eigen::Matrix3d basis;
    Eigen::Vector3d correct;

    // case where the answer can be computed by hand
    double degree = M_PI/180;
    basis << std::cos(80*degree), std::cos(-20*degree), 0.1,
              std::sin(80*degree), std::sin(-20*degree), -100.,
                               0.,                   0.,  0.05;
    Eigen::Matrix3d orth {hexed::custom_math::orthonormal(basis, 2)};
    correct << std::cos(75*degree), std::sin(75*degree), 0.;
    CHECK((orth.col(0) - correct).norm() == Approx(0.).scale(1.));
    correct << std::cos(-15*degree), std::sin(-15*degree), 0.;
    CHECK((orth.col(1) - correct).norm() == Approx(0.).scale(1.));
    correct << 0., 0., 1.;
    CHECK((orth.col(2) - correct).norm() == Approx(0.).scale(1.));

    // test that if you mess with a unitary matrix in the right way, you get the original back
    Eigen::FullPivHouseholderQR<Eigen::Matrix3d> qr (Eigen::Matrix3d::Random());
    orth = qr.matrixQ();
    basis = orth;
    Eigen::Vector3d sum = orth.col(1) + orth.col(2);
    basis.col(1) += 0.2*sum;
    basis.col(2) += 0.2*sum;
    basis.col(0) += Eigen::Vector3d {0.1, -.3, 0.3};
    REQUIRE((hexed::custom_math::orthonormal(basis, 0) - orth).norm() == Approx(0.).scale(1.));

    basis = orth;
    sum = orth.col(0) + orth.col(2);
    basis.col(0) += 3.*sum;
    basis.col(2) += 3.*sum;
    basis.col(1) *= 0.1;
    REQUIRE((hexed::custom_math::orthonormal(basis, 1) - orth).norm() == Approx(0.).scale(1.));
  }
}
