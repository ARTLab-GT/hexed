#include <catch2/catch_all.hpp>
#include <hexed/Occt_geom.hpp>

TEST_CASE("Occt_geom")
{
  REQUIRE_THROWS(hexed::read_cad<IGESControl_Reader>("nonexistent.igs"));
  hexed::Occt_geom geom(hexed::read_cad<IGESControl_Reader>("ellipsoid.igs"));
}
