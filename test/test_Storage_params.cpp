#include <catch2/catch.hpp>

#include <hexed/Storage_params.hpp>

TEST_CASE("Storage_params")
{
  SECTION("3D")
  {
    hexed::Storage_params params {7, 2, 3, 4};
    REQUIRE(params.n_stage == 7);
    REQUIRE(params.n_var == 2);
    REQUIRE(params.n_dim == 3);
    REQUIRE(params.row_size == 4);
    REQUIRE(params.n_qpoint() == 64);
    REQUIRE(params.n_dof() == 128);
    REQUIRE(params.size() == 128*7);
    REQUIRE(params.n_vertices() == 8);
  }
  SECTION("2D")
  {
    hexed::Storage_params params {2, 4, 2, 5};
    REQUIRE(params.size() == 200);
    REQUIRE(params.n_vertices() == 4);
  }
  SECTION("1D")
  {
    hexed::Storage_params params {3, 3, 1, 3};
    REQUIRE(params.size() == 27);
    REQUIRE(params.n_vertices() == 2);
  }
}
