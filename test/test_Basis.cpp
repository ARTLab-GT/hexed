#include "catch.hpp"
#include "Equidistant.hpp"

void test_diff_mat(Basis& basis)
{
  for (int i_result = 0; i_result < basis.rank; ++i_result)
  {
    double derivative = 0;
    for (int i_operand = 0; i_operand < basis.rank; ++i_operand)
    {
      derivative += basis.diff_mat(i_result, i_operand);
    }
    REQUIRE( derivative == Approx(0) );
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
      derivative_lin  += basis.diff_mat(i_result, i_operand)*linear   [i_operand];
      derivative_quad += basis.diff_mat(i_result, i_operand)*quadratic[i_operand];
    }
    REQUIRE( derivative_lin  == Approx(9.07) );
    REQUIRE( derivative_quad == Approx(-0.38 + 2*4.43*basis.node(i_result)) );
  }

  
}

TEST_CASE("Equidistant Basis")
{
  Equidistant equi (4);
  test_diff_mat(equi);
}
