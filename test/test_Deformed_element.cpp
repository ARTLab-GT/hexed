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

  element.jacobian()[0*16] = 1.;
  element.jacobian()[1*16] = 2.;
  element.jacobian()[2*16] = .1;
  element.jacobian()[3*16] = 3.;
  REQUIRE(element.jacobian_determinant(0) == Approx(2.8));

  SECTION("vertex arrangement")
  {
    cartdg::Storage_params params {3, 5, 3, 4};
    cartdg::Deformed_element elem3d {params, {1, 2, -1}, 0.2};
    REQUIRE(elem3d.vertex(0).pos == std::array<double, 3>{0.2, 0.4, -0.2});
    REQUIRE(elem3d.vertex(1).pos == std::array<double, 3>{0.2, 0.4, 0.});
    REQUIRE(elem3d.vertex(2).pos == std::array<double, 3>{0.2, 0.6, -0.2});
    REQUIRE(elem3d.vertex(7).pos == std::array<double, 3>{0.4, 0.6, 0.});
  }
}
