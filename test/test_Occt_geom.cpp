#include <catch2/catch_all.hpp>
#include <hexed/Occt_geom.hpp>

TEST_CASE("Occt_geom")
{
  REQUIRE_THROWS(hexed::Occt_geom::read<IGESControl_Reader>("nonexistent.igs"));
  auto geom = hexed::Occt_geom::read<IGESControl_Reader>("ellipsoid.igs");
  auto nearest = geom.nearest_point(hexed::Mat<3>{-1., 0., 0.});
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.25, 0., 0.}, hexed::math::Approx_equal(0, 1e-12)));
  nearest = geom.nearest_point(hexed::Mat<3>{0., 0., 1.});
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .125}, hexed::math::Approx_equal(0, 1e-12)));
}
