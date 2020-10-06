#include <catch.hpp>

#include <cartdgConfig.hpp>
#include <Equidistant.hpp>
#include <Gauss_lobatto.hpp>

void test_diff_mat(cartdg::Basis& basis)
{
  Eigen::MatrixXd diff_mat = basis.diff_mat();
  for (int i_result = 0; i_result < basis.rank; ++i_result)
  {
    double derivative = 0;
    for (int i_operand = 0; i_operand < basis.rank; ++i_operand)
    {
      derivative += diff_mat(i_result, i_operand);
    }
    REQUIRE( derivative == Approx(0).margin(1e-13) );
  }

  std::vector<double> linear;
  std::vector<double> quadratic;
  for (int i_node = 0; i_node < basis.rank; ++i_node)
  {
    double node = basis.node(i_node);
    linear.push_back(-2.14 + 9.07*node);
    quadratic.push_back(0.07 - 0.38*node - 4.43*node*node);
  }

  for (int i_result = 0; i_result < basis.rank; ++i_result)
  {
    double derivative_lin = 0;
    double derivative_quad = 0;
    for (int i_operand = 0; i_operand < basis.rank; ++i_operand)
    {
      derivative_lin  += diff_mat(i_result, i_operand)*linear   [i_operand];
      derivative_quad += diff_mat(i_result, i_operand)*quadratic[i_operand];
    }
    if (basis.rank > 1)
    {
      REQUIRE( derivative_lin  == Approx(9.07) );
    }
    if (basis.rank > 2)
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
  if (basis.rank >= 3)
  {
    for (int i = 0; i < basis.rank; ++i)
    {
      total += weights(i)*(basis.node(i)*basis.node(i));
    }
    REQUIRE(total == Approx(1./3.));
  }
}

TEST_CASE("Equidistant Basis")
{
  for (int rank = 0; rank < 10; ++rank)
  {
    cartdg::Equidistant equi (rank);
    test_diff_mat(equi);
  }
}

TEST_CASE("Gauss_lobatto Basis")
{
  for (int rank = 2; rank <= MAX_BASIS_RANK; ++rank)
  {
    cartdg::Gauss_lobatto GLo (rank);
    test_diff_mat(GLo);
    test_quadrature(GLo);
  }
}
