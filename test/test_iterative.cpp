#include <catch2/catch_all.hpp>
#include <hexed/iterative.hpp>

TEST_CASE("iterative solvers")
{
  srand(406);
  int n = 10;
  Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(n, n) + 0.1*Eigen::MatrixXd::Random(n, n);
  Eigen::VectorXd rhs = Eigen::VectorXd::Random(n);
  Eigen::VectorXd soln = mat.fullPivLu().solve(rhs);
  hexed::Dense_equation equation(mat, rhs, n*10, Eigen::VectorXd::Random(n));
  SECTION("gmres")
  {
    SECTION("exact") {
      hexed::iterative::gmres(equation, n, 1);
      REQUIRE((equation.vec(0) - soln).norm() == Catch::Approx(0.).scale(1.));
    }
    SECTION("approx") {
      hexed::iterative::gmres(equation, n/2, 2);
      // will probably need to increase the tolerance but i'm curious by how much
      REQUIRE((equation.vec(0) - soln).norm() == Catch::Approx(0.).margin(1e-6));
    }
  }
  SECTION("bicgstab")
  {
    SECTION("exact") {
      hexed::iterative::bicgstab(equation, n);
      REQUIRE((equation.vec(0) - soln).norm() == Catch::Approx(0.).scale(1.));
    }
    SECTION("approx") {
      hexed::iterative::bicgstab(equation, n/2);
      // will probably need to increase the tolerance but i'm curious by how much
      REQUIRE((equation.vec(0) - soln).norm() == Catch::Approx(0.).margin(1e-6));
    }
  }
}
