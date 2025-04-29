#include <catch2/catch_all.hpp>
#include <hexed/brep.hpp>
#include <hexed/Stopwatch.hpp>

TEST_CASE("Line_segment") {
  hexed::Mat<3, 2> endpoints;
  endpoints << 0, 1,
               0, 1,
               1, 1;
  hexed::brep::Line_segment seg(endpoints);
  REQUIRE_THAT(seg.point(hexed::Mat<1>{.1}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0.1, 0.1, 1.}, hexed::math::Approx_equal(0., 1e-12)));
}

TEST_CASE("Circular_arc") {
  hexed::brep::Circular_arc arc({.1, .1, .1}, 10., hexed::constants::pi/2, hexed::constants::pi);
  REQUIRE_THAT(
    arc.point(hexed::Mat<1>{.5}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{.1 - 10*std::sqrt(.5), .1 + 10*std::sqrt(.5), .1},
    hexed::math::Approx_equal())
  );
}

TEST_CASE("Plane") {
  hexed::Mat<3> origin {.1, .1, .1};
  hexed::Mat<3, 2> coords;
  coords <<
    1., 0.,
    0., 1.,
    0., 1.;
  hexed::brep::Plane plane(origin, coords);
  REQUIRE_THAT(plane.point(hexed::Mat<2>{.1, .2}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{.2, .3, .3}, hexed::math::Approx_equal()));
  hexed::Mat<2, 2> bounds;
  bounds << .4, .8,
            .4, 1.;
  plane.reparameterize(bounds);
  REQUIRE_THAT(plane.point(hexed::Mat<2>{.5, .5}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{.7, .8, .8}, hexed::math::Approx_equal()));
}

TEST_CASE("Revolution_surface") {
  hexed::Mat<3, 2> axis_endpoints;
  axis_endpoints <<
    0.01, 0.01,
    0.01, 0.01,
    0.01, 1.01;
  hexed::Mat<3, 2> generatrix_endpoints;
  generatrix_endpoints <<
    1.01, 1.01,
    0.01, 0.01,
    0.01, 1.01;
  hexed::brep::Revolution_surface surf(new hexed::brep::Line_segment(generatrix_endpoints),
                                       hexed::brep::Line_segment(axis_endpoints), 0., hexed::constants::pi);
  REQUIRE_THAT(surf.rotate(hexed::Mat<3>{.01 + std::sqrt(.5), .01 + std::sqrt(.5), .9}, 1.25*hexed::constants::pi),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0.01, -.99, .9}, hexed::math::Approx_equal()));
  REQUIRE_THAT(
    surf.point(hexed::Mat<2>{.5, .5}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{0.01, 1.01, .51}, hexed::math::Approx_equal(0, 1e-3))
  );
}

TEST_CASE("Trimmed_surface") {
  auto plane = std::make_unique<hexed::brep::Plane>(hexed::Mat<3>{4., 4., .0}, hexed::Mat<3, 2>::Identity());
  hexed::Mat<3, 2> endpoints;
  std::vector<hexed::brep::Composite_curve> curves;
  curves.emplace_back();
  endpoints << 0., 3.,
               0., 0.,
               0., 0.;
  curves.back().push_back(std::make_unique<hexed::brep::Line_segment>(endpoints));
  endpoints << 3., 0.,
               0., 3.,
               0., 0.;
  curves.back().push_back(std::make_unique<hexed::brep::Line_segment>(endpoints));
  endpoints << 0., 0.,
               3., 0.,
               0., 0.;
  curves.back().push_back(std::make_unique<hexed::brep::Line_segment>(endpoints));
  curves.emplace_back();
  endpoints << 1., 1.5,
               1., 1.,
               0., 0.;
  curves.back().push_back(std::make_unique<hexed::brep::Line_segment>(endpoints));
  endpoints << 1.5, 1.,
               1., 1.5,
               0., 0.;
  curves.back().push_back(std::make_unique<hexed::brep::Line_segment>(endpoints));
  endpoints << 1., 1.,
               1.5, 1.,
               0., 0.;
  curves.back().push_back(std::make_unique<hexed::brep::Line_segment>(endpoints));
  std::vector<bool> model_space(curves.size(), true);
  hexed::brep::Trimmed_surface trim(plane.release(), std::move(curves), std::move(model_space), 256, 1024);
  // test reparameterization
  REQUIRE_THAT(trim.surface().point(hexed::Mat<2>{0., 0.}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .0}, hexed::math::Approx_equal(0, 1e-6)));
  REQUIRE_THAT(trim.surface().point(hexed::Mat<2>{1., 1.}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{3., 3., .0}, hexed::math::Approx_equal(0, 1e-6)));
  // test is_inside
  CHECK( trim.is_inside(hexed::Mat<2>{.10, .10}));
  CHECK(!trim.is_inside(hexed::Mat<2>{.35, .35}));
  CHECK( trim.is_inside(hexed::Mat<2>{.49, .49}));
  CHECK( trim.is_inside(hexed::Mat<2>{.50, .10}));
  CHECK( trim.is_inside(hexed::Mat<2>{.10, .50}));
  CHECK(!trim.is_inside(hexed::Mat<2>{.55, .55}));
  CHECK(!trim.is_inside(hexed::Mat<2>{-.1, .50}));
  CHECK(!trim.is_inside(hexed::Mat<2>{.50, -.1}));

  // test nearest_point
  REQUIRE_THAT(trim.nearest_point(hexed::Mat<3>{.1, .1, .1}, .2).point(),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, .1, .0}, hexed::math::Approx_equal(0, 1e-3)));
  REQUIRE_THAT(trim.nearest_point(hexed::Mat<3>{1.1, 1.2, .1}, 1.).point(),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{1., 1.2, .0}, hexed::math::Approx_equal(0, 1e-3)));
  REQUIRE_THAT(trim.nearest_point(hexed::Mat<3>{-.1, -.1, .1}, .2).point(),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .0}, hexed::math::Approx_equal(0, 1e-3)));
  REQUIRE(trim.nearest_point(hexed::Mat<3>{.1, .1, .1}, .01).empty());

  // test trimming curves
  const auto& trim_curve = trim.trimming_curves()[1];
  REQUIRE(trim_curve.parameters(0)[0] == Catch::Approx(1.).margin(1e-4));
  REQUIRE(trim_curve.parameters(0)[1] == Catch::Approx(0.).margin(1e-4));
  REQUIRE(trim_curve.parameters(512)[0] == Catch::Approx(.5).margin(1e-4));
  REQUIRE(trim_curve.parameters(512)[1] == Catch::Approx(.5).margin(1e-4));
  REQUIRE(trim_curve.parameters(1024)[0] == Catch::Approx(0.).margin(1e-4));
  REQUIRE(trim_curve.parameters(1024)[1] == Catch::Approx(1.).margin(1e-4));
  CHECK(trim_curve.tangents(0)[0] == Catch::Approx(-std::sqrt(.5)).margin(1e-4));
  CHECK(trim_curve.tangents(0)[1] == Catch::Approx(-std::sqrt(.5)).margin(1e-4));
  CHECK(trim_curve.tangents(0)[2] == Catch::Approx(0.).margin(1e-4));
  CHECK(trim_curve.tangents(512)[0] == Catch::Approx(-std::sqrt(.5)).margin(1e-4));
  CHECK(trim_curve.tangents(512)[1] == Catch::Approx(-std::sqrt(.5)).margin(1e-4));
  CHECK(trim_curve.tangents(512)[2] == Catch::Approx(0.).margin(1e-4));
  CHECK(trim_curve.tangents(1023)[0] == Catch::Approx(-std::sqrt(.5)).margin(1e-4));
  CHECK(trim_curve.tangents(1023)[1] == Catch::Approx(-std::sqrt(.5)).margin(1e-4));
  CHECK(trim_curve.tangents(1023)[2] == Catch::Approx(0.).margin(1e-4));
}

TEST_CASE("Geom_3d", "[.slow]") {
  #ifdef DEBUG
  hexed::Int n_div_min = 8;
  hexed::Int n_div_max = 128;
  bool vis_volume = false;
  #else
  hexed::Int n_div_min = 128;
  hexed::Int n_div_max = 1024;
  bool vis_volume = true;
  #endif
  SECTION("cylinder_extruded") {
    hexed::brep::Geom_3d geom("../test_assets/cylinder_extruded.iges", n_div_min, n_div_max);
    geom.visualize("default", "cylinder_extruded", 100, vis_volume);
  }
  #if 0
  SECTION("prism_twisted") {
    hexed::brep::Geom_3d geom("../test_assets/prism_twisted.iges", n_div_min, n_div_max);
    geom.visualize("default", "prism_twisted", 30, vis_volume);
  }
  SECTION("weird_surface") {
    hexed::brep::Geom_3d geom("../test_assets/weird_surface.iges", n_div_min, n_div_max);
    hexed::Mat<3, 2> bounds;
    bounds <<
      -.1, .1,
      -.1, .1,
       0., .2;
    geom.visualize("default", "weird_surface", 30, vis_volume, bounds);
  }
  SECTION("misleading_normal") {
    hexed::Stopwatch sw;
    sw.start();
    hexed::brep::Geom_3d geom("../test_assets/misleading_normal.iges", n_div_min, n_div_max);
    sw.pause();
    std::cout << "startup time: " << sw.time() << std::endl;
    sw.reset();
    sw.start();
    auto& surf = geom.surfaces()[0];
    hexed::Int n = hexed::math::pow(10, 5);
    #pragma omp parallel for
    for (hexed::Int i = 0; i < n; ++i) {
      // we're going to do something with the results just to make sure the calculation isn't optimized away
      HEXED_ASSERT(std::isfinite(surf.point((hexed::Mat<2>::Random() + hexed::Mat<2>::Ones())/2).squaredNorm()), "")
    }
    sw.pause();
    std::cout << "point evaluation time: " << sw.time()/n << "s/point" << std::endl;
    sw.reset();
    hexed::Mat<3, 2> bounds;
    bounds <<
      -.2, .2,
      -.1, .3,
      -.2, .2;
    sw.start();
    geom.visualize("default", "misleading_normal", 30, vis_volume, bounds);
    sw.pause();
    std::cout << "visualization time: " << sw.time() << std::endl;
    sw.reset();
  }
  #endif
}
