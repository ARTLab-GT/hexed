#include <catch2/catch_all.hpp>
#include <hexed/cad_geom.hpp>

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
  hexed::cad_geom::Revolution_surface surf(new hexed::cad_geom::Line_segment(generatrix_endpoints),
                                           hexed::cad_geom::Line_segment(axis_endpoints), 0., hexed::constants::pi);
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
  REQUIRE(n(2) == Catch::Approx(.1));
}

TEST_CASE("Geom") {
  hexed::cad_geom::Geom geom("../test_assets/cylinder_extruded.iges");
  geom.visualize("cylinder_extruded");
}
