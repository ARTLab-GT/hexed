#include <catch2/catch_all.hpp>
#include <hexed/constants.hpp>
#include <hexed/standard_atmosphere.hpp>

TEST_CASE("standard_atmosphere")
{
  auto state = hexed::standard_atmosphere_geopot(0);
  REQUIRE(state[0] == Catch::Approx(1.225));
  REQUIRE(state[1] == Catch::Approx(hexed::constants::atmosphere));
  state = hexed::standard_atmosphere_geopot(1e4, -10);
  REQUIRE(state[0] == Catch::Approx(0.4320687958));
  REQUIRE(state[1] == Catch::Approx(26436.26759));
  state = hexed::standard_atmosphere(1e4, -10);
  REQUIRE(state[0] == Catch::Approx(0.4320687958).epsilon(1e-2)); // geometric altitude version should be approximately the same as geopotential
  REQUIRE(state[1] == Catch::Approx(26436.26759).epsilon(1e-2));
  state = hexed::standard_atmosphere_geopot(8e4, 0);
  REQUIRE(state[0] == Catch::Approx(0.00001570054961));
  REQUIRE(state[1] == Catch::Approx(0.8862795041));
}
