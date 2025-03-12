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
  hexed::Mat<3, 2> points;
  points << 0., .5,
            .6, .6,
            1., 1.;
  auto inter_points = seg.intersection_params(points);
  REQUIRE(inter_points.size() == 1);
  REQUIRE(inter_points[0].params(0) == Catch::Approx(.6));
  REQUIRE(inter_points[0].interp_coef == Catch::Approx(1.2));
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
  hexed::Mat<3, 2> points;
  points <<  -4.9,  -4.9,
             10.1,    .1,
               .1,    .1;
  auto inter_points = arc.intersection_params(points);
  REQUIRE(inter_points.size() == 1);
  REQUIRE(inter_points[0].params(0) == Catch::Approx(1./3.));
  REQUIRE(inter_points[0].interp_coef == Catch::Approx(1. - std::sqrt(.75)));
  points << -10.9,    .1,
               .1,  11.1,
               .1,    .1;
  REQUIRE(arc.intersection_params(points).size() == 2);
  points << -19.9,    .1,
               .1,  20.1,
               .1,    .1;
  REQUIRE(arc.intersection_params(points).size() == 0);
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
  REQUIRE_THAT(plane.nearest_point(hexed::Mat<3>{.5, 0.1, 1.1}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{.5, .6, .6}, hexed::math::Approx_equal()));
  hexed::Mat<3, 2> endpoints;
  endpoints <<
    1., 1.,
    .5, .5,
    2., 0.;
  auto sects = plane.intersection_params(endpoints);
  REQUIRE(sects.size() == 1);
  REQUIRE(sects[0].params[0] == Catch::Approx(.9));
  REQUIRE(sects[0].params[1] == Catch::Approx(.4));
  REQUIRE(sects[0].interp_coef == Catch::Approx(.75));
  endpoints <<
    -1., -1.,
    .5, .5,
    2., 0.;
  REQUIRE(plane.intersection_params(endpoints).size() == 0);
  endpoints <<
    1., 1.,
    5., 5.,
    2., 0.;
  REQUIRE(plane.intersection_params(endpoints).size() == 0);
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
                                       hexed::brep::Line_segment(axis_endpoints), 1024, 0., hexed::constants::pi);
  REQUIRE_THAT(surf.rotate(hexed::Mat<3>{.01 + std::sqrt(.5), .01 + std::sqrt(.5), .9}, 1.25*hexed::constants::pi),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0.01, -.99, .9}, hexed::math::Approx_equal()));
  REQUIRE_THAT(
    surf.point(hexed::Mat<2>{.5, .5}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{0.01, 1.01, .51}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{0.01, 1.01, .51}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{0.01, 1.01, .51}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{1.01, 1.01, .71}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{std::sqrt(.5) + .01, std::sqrt(.5) + .01, .71},
                                 hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{-.99, -.99, .71}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{-.99, 0.01, .71}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{0.01, 1.01, 2.01}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{0.01, 1.01, 1.01}, hexed::math::Approx_equal(0, 1e-3))
  );
  REQUIRE_THAT(
    surf.nearest_point(hexed::Mat<3>{1.01, -.99, -.11}),
    Catch::Matchers::RangeEquals(hexed::Mat<3>{1.01, 0.01, 0.01}, hexed::math::Approx_equal(0, 1e-3))
  );
  hexed::Mat<3> n = surf.nearest_point(hexed::Mat<3>{0.01, 0.01, .11});
  REQUIRE(!std::isnan(n(0)));
  REQUIRE(!std::isnan(n(1)));
  REQUIRE(n(2) == Catch::Approx(.11).epsilon(1e-2));
  SECTION("intersections") { // note: the following calculations should be exact up to rounding errors
    hexed::Mat<3, 2> endpoints;
    endpoints <<
      -0.99, 1.01,
      .51, .51,
      .01, 1.01;
    auto sections = surf.intersection_params(endpoints);
    REQUIRE(sections.size() == 2);
    std::vector<double> test {sections[0].params(0), sections[1].params(0)};
    std::vector<double> correct {.5*(1 - std::sqrt(.75)), .5*(1 + std::sqrt(.75))};
    REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
    test = std::vector<double>{sections[0].interp_coef, sections[1].interp_coef};
    correct = std::vector<double>{.5*(1 - std::sqrt(.75)), .5*(1 + std::sqrt(.75))};
    REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
    test = std::vector<double>{sections[0].params(1), sections[1].params(1)};
    correct = std::vector<double>{5./6., 1./6.};
    REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
    SECTION("horizontal section line and only 1 intersection") {
      endpoints <<
        .51, .51,
        -.99, 1.01,
        1./7. + .01, 1./7. + .01;
      auto sections = surf.intersection_params(endpoints);
      REQUIRE(sections.size() == 1);
      std::vector<double> test {sections[0].params(0)};
      std::vector<double> correct {1./7.};
      REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
      test = std::vector<double>{sections[0].interp_coef};
      correct = std::vector<double>{.5*(1 + std::sqrt(.75))};
      REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
    }
    SECTION("annular surface") {
      hexed::Mat<3, 2> new_generatrix_endpoints;
      new_generatrix_endpoints <<
        1.01, 2.01,
        0.01, 0.01,
        0.01, 0.01;
      hexed::brep::Revolution_surface new_surf(new hexed::brep::Line_segment(new_generatrix_endpoints),
                                               hexed::brep::Line_segment(axis_endpoints), 1024,
                                               0., hexed::constants::pi);
      endpoints <<
        1.01, 1.01,
        1.01, 1.01,
        1.00, 0.00;
      auto sections = new_surf.intersection_params(endpoints);
      REQUIRE(sections.size() == 1);
      std::vector<double> test {sections[0].params(0)};
      std::vector<double> correct {std::sqrt(2.) - 1.};
      REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
      test = std::vector<double>{sections[0].interp_coef};
      correct = std::vector<double>{.99};
      REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
      test = std::vector<double>{sections[0].params(1)};
      correct = std::vector<double>{.25};
      REQUIRE_THAT(test, Catch::Matchers::UnorderedRangeEquals(correct, hexed::math::Approx_equal(0, 1e-6)));
    }
  }
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
  hexed::brep::Trimmed_surface trim(plane.release(), std::move(curves), std::move(model_space), 1024);
  // test reparameterization
  REQUIRE_THAT(trim.surface().point(hexed::Mat<2>{0., 0.}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .0}, hexed::math::Approx_equal(0, 1e-6)));
  REQUIRE_THAT(trim.surface().point(hexed::Mat<2>{1., 1.}),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{3., 3., .0}, hexed::math::Approx_equal(0, 1e-6)));
  // test is_inside
  REQUIRE( trim.is_inside(hexed::Mat<2>{.10, .10}));
  REQUIRE(!trim.is_inside(hexed::Mat<2>{.35, .35}));
  REQUIRE( trim.is_inside(hexed::Mat<2>{.49, .49}));
  REQUIRE( trim.is_inside(hexed::Mat<2>{.50, .10}));
  REQUIRE( trim.is_inside(hexed::Mat<2>{.10, .50}));
  REQUIRE(!trim.is_inside(hexed::Mat<2>{.55, .55}));
  REQUIRE(!trim.is_inside(hexed::Mat<2>{-.1, .50}));
  REQUIRE(!trim.is_inside(hexed::Mat<2>{.50, -.1}));

  // test nearest_point
  REQUIRE_THAT(trim.nearest_point(hexed::Mat<3>{.1, .1, .1}, .2).point(),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, .1, .0}, hexed::math::Approx_equal(0, 1e-3)));
  REQUIRE_THAT(trim.nearest_point(hexed::Mat<3>{1.1, 1.2, .1}, 1.).point(),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{1., 1.2, .0}, hexed::math::Approx_equal(0, 1e-3)));
  REQUIRE_THAT(trim.nearest_point(hexed::Mat<3>{-.1, -.1, .1}, .2).point(),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .0}, hexed::math::Approx_equal(0, 1e-3)));
  REQUIRE(trim.nearest_point(hexed::Mat<3>{.1, .1, .1}, .01).empty());
}

TEST_CASE("Geom_3d", "[.slow]") {
  #ifdef DEBUG
  hexed::Int n_div = 128;
  bool vis_volume = false;
  #else
  hexed::Int n_div = 1024;
  bool vis_volume = true;
  #endif
  #if 0
  SECTION("cylinder_extruded") {
    hexed::brep::Geom_3d geom("../test_assets/cylinder_extruded.iges", n_div);
    geom.visualize("default", "cylinder_extruded", 100, vis_volume);
  }
  #endif
  SECTION("prism_twisted") {
    hexed::brep::Geom_3d geom("../test_assets/prism_twisted.iges", n_div);
    geom.visualize("default", "prism_twisted", 30, vis_volume);
  }
  #if 0
  SECTION("weird_surface") {
    hexed::brep::Geom_3d geom("../test_assets/weird_surface.iges", n_div);
    hexed::Mat<3, 2> bounds;
    bounds <<
      -.1, .1,
      -.1, .1,
      0., .2;
    geom.visualize("default", "weird_surface", 30, vis_volume, bounds);
  }
  #endif
}
