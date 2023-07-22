#include <catch2/catch_all.hpp>
#include <hexed/Occt_geom.hpp>
#include <hexed/constants.hpp>
#include <hexed/Simplex_geom.hpp>

#if HEXED_USE_OCCT

void test(std::string file_extension)
{
  REQUIRE_THROWS(hexed::Occt_geom::read("nonexistent." + file_extension));
  hexed::Occt_geom geom(hexed::Occt_geom::read("ellipsoid." + file_extension), 3);
  auto nearest = geom.nearest_point(-hexed::Mat<3>::Unit(0));
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.25, 0., 0.}, hexed::math::Approx_equal(0, 1e-12)));
  nearest = geom.nearest_point(hexed::Mat<3>::Unit(2));
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .125}, hexed::math::Approx_equal(0, 1e-12)));
  auto intersections = geom.intersections(hexed::Mat<3>::Unit(2), 3*hexed::Mat<3>::Unit(2));
  REQUIRE_THAT(intersections, Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-.4375, -.5625}, hexed::math::Approx_equal(0, 1e-12)));
}

TEST_CASE("Occt_geom")
{
  SECTION("3D")
  {
    test("igs");
    test("stp");
  }
  SECTION("2D")
  {
    auto shape = hexed::Occt_geom::read("ellipse.STEP");
    hexed::Occt_geom::write_image(shape, "test.png");
    std::vector<std::unique_ptr<hexed::Surface_geom>> geoms;
    geoms.emplace_back(new hexed::Occt_geom(shape, 2));
    geoms.emplace_back(new hexed::Simplex_geom(hexed::segments(shape, 1000)));
    for (auto& geom : geoms) {
      auto nearest = geom->nearest_point(-hexed::Mat<2>::Unit(0));
      REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<2>{-.25, 0.}, hexed::math::Approx_equal(0, 1e-12)));
      nearest = geom->nearest_point(hexed::Mat<2>::Unit(1));
      REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<2>{0., .125}, hexed::math::Approx_equal(0, 1e-12)));
      auto intersections = geom->intersections(hexed::Mat<2>::Unit(1), 3*hexed::Mat<2>::Unit(1));
      REQUIRE_THAT(intersections, Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-.4375, -.5625}, hexed::math::Approx_equal(0, 1e-12)));
    }
  }
}

#endif
