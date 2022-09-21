#include <catch2/catch.hpp>
#include <hexed/config.hpp>
#include <hexed/Derivative.hpp>
#include <hexed/Gauss_legendre.hpp>

double arbitrary_polynomial(double pos, int i)
{
  return 0.1*pos*pos - i*pos - 3.;
}

TEST_CASE("Derivative")
{
  const int row_size = hexed::config::max_row_size;
  static_assert (row_size >= 3); // this test requires at least quadratic basis
  hexed::Gauss_legendre basis (row_size);
  Eigen::Matrix<double, row_size, 3> qpoint_vals;
  Eigen::Matrix<double, 2, 3> boundary_vals;
  for (int i_var = 0; i_var < 3; ++i_var) {
    // write some arbitrary quadratic polynomial
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      qpoint_vals(i_qpoint, i_var) = arbitrary_polynomial(basis.node(i_qpoint), i_var);
    }
    // set the boundary values to match the interior
    for (int i_side : {0, 1}) {
      boundary_vals(i_side, i_var) = arbitrary_polynomial(i_side, i_var);
    }
  }
  // mess with the boundary values of the last variable so we can check conservation
  boundary_vals(0, 2) = 0.2;
  boundary_vals(1, 2) = 0.5;
  hexed::Derivative<row_size> derivative (basis);
  auto result = derivative(qpoint_vals, boundary_vals);
  // test that the derivative is correct when the boundary values match the interior
  for (int i_var = 0; i_var < 2; ++i_var) {
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      REQUIRE(result(i_qpoint, i_var) == Approx(0.2*basis.node(i_qpoint) - i_var).scale(1.));
    }
  }
  // check that the integral of the derivative is the difference between the boundary values,
  // even when this disagrees with the derivative of the interior (i.e. the derivative is conservative)
  double integral = 0.;
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    integral += result(i_qpoint, 2)*basis.node_weights()(i_qpoint);
  }
  REQUIRE(integral == Approx(0.5 - 0.2));
}
