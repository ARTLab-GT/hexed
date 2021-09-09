#include <catch2/catch.hpp>

#include <Deformed_element.hpp>

TEST_CASE("Deformed_element.hpp")
{
  cartdg::Storage_params params {2, 2, 2, 4};
  cartdg::Deformed_element element {params};
  for (int i_jac = 0; i_jac < 2*2*16; ++i_jac)
  {
    element.jacobian()[i_jac] = i_jac/4.;
  }
  for (int i_dim = 0, i_jac = 0; i_dim < 2; ++i_dim)
  {
    for (int j_dim = 0; j_dim < 2; ++j_dim)
    {
      for (int i_qpoint = 0; i_qpoint < 16; ++i_qpoint)
      {
        REQUIRE(element.jacobian(i_dim, j_dim, i_qpoint) == i_jac++/4.);
      }
    }
  }
}
