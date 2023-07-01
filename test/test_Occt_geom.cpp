#include <catch2/catch_all.hpp>
#include <hexed/Occt_geom.hpp>
#include <hexed/constants.hpp>

TEST_CASE("Occt_geom")
{
  REQUIRE_THROWS(hexed::Occt_geom::read<IGESControl_Reader>("nonexistent.igs"));
  auto geom = hexed::Occt_geom::read<IGESControl_Reader>("ellipsoid.igs");
  auto nearest = geom.nearest_point(-hexed::Mat<3>::Unit(0));
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.25, 0., 0.}, hexed::math::Approx_equal(0, 1e-12)));
  nearest = geom.nearest_point(hexed::Mat<3>::Unit(2));
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .125}, hexed::math::Approx_equal(0, 1e-12)));
  auto intersections = geom.intersections(hexed::Mat<3>::Unit(2), 3*hexed::Mat<3>::Unit(2));
  REQUIRE_THAT(intersections, Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-.4375, -.5625}, hexed::math::Approx_equal(0, 1e-12)));
}
