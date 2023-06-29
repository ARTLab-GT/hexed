#include <catch2/catch_all.hpp>
#include <hexed/Occt_geom.hpp>

TEST_CASE("Occt_geom")
{
  REQUIRE_THROWS(hexed::Occt_geom::read<IGESControl_Reader>("nonexistent.igs"));
  auto geom = hexed::Occt_geom::read<IGESControl_Reader>("ellipsoid.igs");
}
