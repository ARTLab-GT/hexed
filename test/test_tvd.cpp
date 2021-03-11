#include <catch.hpp>

#include <cartdgConfig.hpp>
#include <kernels/local/n_extrema.hpp>

TEST_CASE("n_extrema")
{
  double read [7] {};
  REQUIRE(cartdg::n_extrema<0>(read) == 0);
  REQUIRE(cartdg::n_extrema<1>(read) == 0);
  read[1] = 1.;
  REQUIRE(cartdg::n_extrema<2>(read) == 0);
  REQUIRE(cartdg::n_extrema<3>(read) == 1);
  read[1] = 0;
  REQUIRE(cartdg::n_extrema<3>(read) == 0);
  REQUIRE(cartdg::n_extrema<6>(read) == 0);
  read[1] = 1;
  read[2] = -1;
  REQUIRE(cartdg::n_extrema<6>(read) == 2);
  read[4] = -1;
  REQUIRE(cartdg::n_extrema<6>(read) == 4);
  read[3] = -2;
  REQUIRE(cartdg::n_extrema<6>(read) == 2);
}