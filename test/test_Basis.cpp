#include <catch2/catch_all.hpp>

#include <hexed/config.hpp>
#include <hexed/Equidistant.hpp>
#include <hexed/Gauss_lobatto.hpp>
#include <hexed/Gauss_legendre.hpp>

void test_diff_mat(hexed::Basis& basis)
{
  Eigen::MatrixXd diff_mat = basis.diff_mat();
  for (int i_result = 0; i_result < basis.row_size; ++i_result)
  {
    double derivative = 0;
    for (int i_operand = 0; i_operand < basis.row_size; ++i_operand)
    {
      derivative += diff_mat(i_result, i_operand);
    }
    REQUIRE( derivative == Catch::Approx(0).margin(1e-13) );
  }

  std::vector<double> linear;
  std::vector<double> quadratic;
  for (int i_node = 0; i_node < basis.row_size; ++i_node)
  {
    double node = basis.node(i_node);
    linear.push_back(-2.14 + 9.07*node);
    quadratic.push_back(0.07 - 0.38*node - 4.43*node*node);
  }

  for (int i_result = 0; i_result < basis.row_size; ++i_result)
  {
    double derivative_lin = 0;
    double derivative_quad = 0;
    for (int i_operand = 0; i_operand < basis.row_size; ++i_operand)
    {
      derivative_lin  += diff_mat(i_result, i_operand)*linear   [i_operand];
      derivative_quad += diff_mat(i_result, i_operand)*quadratic[i_operand];
    }
    if (basis.row_size > 1)
    {
      REQUIRE( derivative_lin  == Catch::Approx(9.07) );
    }
    if (basis.row_size > 2)
    {
      REQUIRE( derivative_quad == Catch::Approx(-0.38 - 2*4.43*basis.node(i_result)) );
    }
  }
}

void test_quadrature(hexed::Basis& basis)
{
  Eigen::VectorXd weights = basis.node_weights();
  REQUIRE(weights.sum() == Catch::Approx(1.));
  double total = 0.;
  if (basis.row_size >= 3)
  {
    for (int i = 0; i < basis.row_size; ++i)
    {
      total += weights(i)*(basis.node(i)*basis.node(i));
    }
    REQUIRE(total == Catch::Approx(1./3.));
  }
}

void test_boundary(hexed::Basis& basis)
{
  auto boundary = basis.boundary();
  REQUIRE(boundary.rows() == 2);
  REQUIRE(boundary.cols() == basis.row_size);
  Eigen::VectorXd node_vals (basis.row_size);
  for (int i_node = 0; i_node < basis.row_size; ++i_node)
  {
    node_vals[i_node] = 0.15*basis.node(i_node) + 0.37;
  }
  Eigen::VectorXd boundary_vals = boundary*node_vals;
  REQUIRE(boundary_vals(0) == Catch::Approx(0.37));
  REQUIRE(boundary_vals(1) == Catch::Approx(0.52));
}

void test_orthogonal(hexed::Basis& basis)
{
  Eigen::VectorXd weights = basis.node_weights();
  for (int i_orth = 0; i_orth < basis.row_size; ++i_orth)
  {
    Eigen::VectorXd weighted_orth = basis.orthogonal(i_orth).cwiseProduct(weights);
    for (int j_orth = 0; j_orth < basis.row_size; ++j_orth)
    {
      auto orth = basis.orthogonal(j_orth);
      REQUIRE(weighted_orth.dot(orth) == Catch::Approx(i_orth == j_orth ? 1. : 0.).margin(1e-10));
    }
  }
}

void test_transform(hexed::Basis& basis)
{
  const int rs {basis.row_size};
  const double tol {1e-10};

  // Prolongation followed by restriction should be identity operator.
  Eigen::MatrixXd id {Eigen::MatrixXd::Zero(rs, rs)};
  for (int i_half : {0, 1}) id += basis.restrict(i_half)*basis.prolong(i_half);
  REQUIRE((Eigen::MatrixXd::Identity(rs, rs) - id).norm() < tol);

  // Restricting and then prolonging monomial basis should yield monomial basis.
  Eigen::MatrixXd prolong  (2*rs, rs);
  Eigen::MatrixXd restrict (rs, 2*rs);
  Eigen::MatrixXd coarse_monomial (2*rs, rs);
  for (int i_half : {0, 1}) {
    prolong.block(i_half*rs, 0, rs, rs) = basis.prolong(i_half);
    restrict.block(0, i_half*rs, rs, rs) = basis.restrict(i_half);
    for (int degree = 0; degree < rs; ++degree) {
      for (int i_node = 0; i_node < rs; ++i_node) {
        coarse_monomial(i_half*rs + i_node, degree) = std::pow((basis.node(i_node) + i_half)/2., degree);
      }
    }
  }
  Eigen::MatrixXd same {prolong*restrict*coarse_monomial};
  CHECK((same - coarse_monomial).norm() < tol);
}

TEST_CASE("Equidistant Basis")
{
  for (int row_size = 2; row_size < 10; ++row_size)
  {
    hexed::Equidistant equi (row_size);
    test_diff_mat(equi);
    test_boundary(equi);
  }
}

TEST_CASE("Gauss_lobatto Basis")
{
  for (int row_size = 2; row_size <= hexed::config::max_row_size; ++row_size)
  {
    hexed::Gauss_lobatto GLo (row_size);
    test_diff_mat(GLo);
    test_quadrature(GLo);
    test_orthogonal(GLo);
    test_boundary(GLo);
  }
}

TEST_CASE("Gauss_legendre Basis")
{
  for (int row_size = 2; row_size <= hexed::config::max_row_size; ++row_size)
  {
    hexed::Gauss_legendre GLe (row_size);
    test_diff_mat(GLe);
    test_quadrature(GLe);
    test_orthogonal(GLe);
    test_boundary(GLe);
    test_transform(GLe);
  }
}

TEST_CASE("interpolation")
{
  hexed::Equidistant basis {5};
  Eigen::VectorXd values {5};
  for (int i = 0; i < 5; ++i) values[i] = std::pow(i/4., 3);
  SECTION("large sample")
  {
    Eigen::VectorXd sample {{0, 0.5, 0.7, 0.16, 1.2, 0., 0.}};
    Eigen::VectorXd interpolated = basis.interpolate(sample)*values;
    REQUIRE(interpolated[0] == 0.);
    REQUIRE(interpolated[1] == Catch::Approx(0.125));
    REQUIRE(interpolated[2] == Catch::Approx(.7*.7*.7));
    REQUIRE(interpolated[4] == Catch::Approx(1.2*1.2*1.2));
    REQUIRE(interpolated[6] == 0.);
  }
  SECTION("small sample")
  {
    Eigen::VectorXd sample {{0.2}};
    Eigen::VectorXd interpolated = basis.interpolate(sample)*values;
    REQUIRE(interpolated[0] == Catch::Approx(.2*.2*.2));
  }
}
