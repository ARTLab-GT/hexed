#include <catch2/catch.hpp>

#include <Storage_params.hpp>

TEST_CASE("Storage_params")
{
  SECTION("3D")
  {
    cartdg::Storage_params params {7, 2, 3, 4};
    REQUIRE(params.n_stage == 7);
    REQUIRE(params.n_var == 2);
    REQUIRE(params.n_dim == 3);
    REQUIRE(params.row_size == 4);
    REQUIRE(params.n_qpoint() == 64);
    REQUIRE(params.n_dof() == 128);
    REQUIRE(params.size() == 128*7);
  }
  SECTION("2D")
  {
    cartdg::Storage_params params {2, 4, 2, 5};
    REQUIRE(params.size() == 200);
  }
  SECTION("1D")
  {
    cartdg::Storage_params params {3, 3, 1, 3};
    REQUIRE(params.size() == 27);
  }
}
