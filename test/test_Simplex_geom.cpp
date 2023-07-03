#include <catch2/catch_all.hpp>
#include <hexed/Simplex_geom.hpp>

TEST_CASE("Simplex_geom")
{
  SECTION("static 2d")
  {
    std::vector<hexed::Mat<2, 2>> elems(2);
    elems[0] << 0, 2,
                1, 1;
    elems[1] << 2, 2,
                0, 1;
    hexed::Simplex_geom<2> geom(elems);
    hexed::Mat<> nearest;
    nearest = geom.nearest_point(hexed::Mat<2>{.5, 2.});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<2>{.5, 1.}, hexed::math::Approx_equal(0, 1e-12)));
    nearest = geom.nearest_point(hexed::Mat<2>{3, .5});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<2>{2., .5}, hexed::math::Approx_equal(0, 1e-12)));
    nearest = geom.nearest_point(hexed::Mat<2>{3, 3});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<2>{2., 1.}, hexed::math::Approx_equal(0, 1e-12)));
  }
  SECTION("dynamic 3d")
  {
    std::vector<hexed::Mat<3, 3>> elems;
    elems.emplace_back(3, 3);
    elems[0] << 1, 1, 0,
                0, 1, 0,
                0, 0, 1;
    hexed::Simplex_geom<3> geom(elems);
    hexed::Mat<> nearest;
    // nearest on face
    nearest = geom.nearest_point(hexed::Mat<3>{1., .1, 1.});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{.5, .1, .5}, hexed::math::Approx_equal(0, 1e-12)));
    // nearest on edges
    nearest = geom.nearest_point(hexed::Mat<3>{2., .1, 0.});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{1., .1, 0.}, hexed::math::Approx_equal(0, 1e-12)));
    nearest = geom.nearest_point(hexed::Mat<3>{0., 1.5, 1.});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{.5, .5, .5}, hexed::math::Approx_equal(0, 1e-12)));
    nearest = geom.nearest_point(hexed::Mat<3>{.5, -.1, .5});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{.5, 0., .5}, hexed::math::Approx_equal(0, 1e-12)));
    // nearest at vertices
    nearest = geom.nearest_point(hexed::Mat<3>{1.1, 2., 0.1});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{1., 1., 0.}, hexed::math::Approx_equal(0, 1e-12)));
    nearest = geom.nearest_point(hexed::Mat<3>{-.1, -.1, 1.1});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., 1.}, hexed::math::Approx_equal(0, 1e-12)));
    nearest = geom.nearest_point(hexed::Mat<3>{1., -.1, -.1});
    REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{1., 0., 0.}, hexed::math::Approx_equal(0, 1e-12)));
  }
}
