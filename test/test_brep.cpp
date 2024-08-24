#include <catch2/catch_all.hpp>
#include <hexed/brep.hpp>

TEST_CASE("Line_segment") {
  hexed::Mat<3, 2> endpoints;
  endpoints << 0, 1,
               0, 1,
               1, 1;
  hexed::brep::Line_segment seg(endpoints);
  REQUIRE_THAT(seg.point(hexed::Mat<1>{.1}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0.1, 0.1, 1.}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE_THAT(seg.nearest_point(hexed::Mat<3>{1., 0., 0.}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0.5, 0.5, 1.}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE_THAT(seg.nearest_point(hexed::Mat<3>{-1., 0., 1.}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., 1.}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE_THAT(seg.nearest_point(hexed::Mat<3>{2., 1., 1.}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{1., 1., 1.}, hexed::math::Approx_equal(0., 1e-12)));
}

TEST_CASE("Circular_arc") {
  hexed::brep::Circular_arc arc({.1, .1, .1}, 10., hexed::constants::pi/2, hexed::constants::pi);
  REQUIRE_THAT(
    arc.point(hexed::Mat<1>{.5}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{.1 - 10*std::sqrt(.5), .1 + 10*std::sqrt(.5), .1},
    hexed::math::Approx_equal())
  );
  REQUIRE_THAT(
    arc.nearest_point(hexed::Mat<3>{-9.9, 10.1, 5.}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{.1 - 10*std::sqrt(.5), .1 + 10*std::sqrt(.5), .1},
    hexed::math::Approx_equal())
  );
  REQUIRE_THAT(
    arc.nearest_point(hexed::Mat<3>{10.1, 10.1, 5.}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, 10.1, .1},
    hexed::math::Approx_equal())
  );
  REQUIRE_THAT(
    arc.nearest_point(hexed::Mat<3>{-9.9, -10.1, 5.}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{-9.9, .1, .1},
    hexed::math::Approx_equal())
  );
}

TEST_CASE("Revolution_surface") {
  hexed::Mat<3, 2> axis_endpoints;
  axis_endpoints <<
    0., 0.,
    0., 0.,
    0., 1.;
  hexed::Mat<3, 2> generatrix_endpoints;
  generatrix_endpoints <<
    1., 1.,
    0., 0.,
    0., 1.;
  hexed::brep::Revolution_surface surf(new hexed::brep::Line_segment(generatrix_endpoints),
                                           hexed::brep::Line_segment(axis_endpoints), 0., hexed::constants::pi);
  REQUIRE_THAT(
    surf.point(hexed::Mat<2>{.5, .5}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 1., .5}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{0., 1., .5}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 1., .5}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{1., 1., .7}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{std::sqrt(.5), std::sqrt(.5), .7}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{-1., -1., .7}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{-1., 0., .7}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{0., 1., 2.}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 1., 1.}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{1., -1., -.1}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{1., 0., 0.}, hexed::math::Approx_equal(0, 1e-3))
  );
  hexed::Mat<3> n = surf.nearest_point(hexed::Mat<3>{0., 0., .1});
  REQUIRE(!std::isnan(n(0)));
  REQUIRE(!std::isnan(n(1)));
  REQUIRE(n(2) == Catch::Approx(.1).epsilon(1e-2));
}

TEST_CASE("Geom", "[.slow]") {
  hexed::brep::Geom geom("../test_assets/cylinder_extruded.iges");
  geom.visualize("cylinder_extruded");
}
