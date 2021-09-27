#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Equidistant.hpp>
#include <Gauss_lobatto.hpp>
#include <Gauss_legendre.hpp>

void test_diff_mat(cartdg::Basis& basis)
{
  Eigen::MatrixXd diff_mat = basis.diff_mat();
  for (int i_result = 0; i_result < basis.row_size; ++i_result)
  {
    double derivative = 0;
    for (int i_operand = 0; i_operand < basis.row_size; ++i_operand)
    {
      derivative += diff_mat(i_result, i_operand);
    }
    REQUIRE( derivative == Approx(0).margin(1e-13) );
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
      REQUIRE( derivative_lin  == Approx(9.07) );
    }
    if (basis.row_size > 2)
    {
      REQUIRE( derivative_quad == Approx(-0.38 - 2*4.43*basis.node(i_result)) );
    }
  }
}

void test_quadrature(cartdg::Basis& basis)
{
  Eigen::VectorXd weights = basis.node_weights();
  REQUIRE(weights.sum() == Approx(1.));
  double total = 0.;
  if (basis.row_size >= 3)
  {
    for (int i = 0; i < basis.row_size; ++i)
    {
      total += weights(i)*(basis.node(i)*basis.node(i));
    }
    REQUIRE(total == Approx(1./3.));
  }
}

TEST_CASE("Equidistant Basis")
{
  for (int row_size = 0; row_size < 10; ++row_size)
  {
    cartdg::Equidistant equi (row_size);
    test_diff_mat(equi);
  }
}

TEST_CASE("Gauss_lobatto Basis")
{
  for (int row_size = 2; row_size <= CARTDG_MAX_BASIS_ROW_SIZE; ++row_size)
  {
    cartdg::Gauss_lobatto GLo (row_size);
    test_diff_mat(GLo);
    test_quadrature(GLo);
  }
}

TEST_CASE("Gauss_legendre Basis")
{
  for (int row_size = 2; row_size <= CARTDG_MAX_BASIS_ROW_SIZE; ++row_size)
  {
    cartdg::Gauss_legendre GLe (row_size);
    test_diff_mat(GLe);
    test_quadrature(GLe);
  }
}
