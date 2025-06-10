#include <catch2/catch_all.hpp>
#include <hexed/Surface_geom.hpp>
#include <hexed/Tree_curve_edge.hpp>

TEST_CASE("Compound_edge") {
  hexed::Array<double> nodes({4, 3});
  nodes(0) = 0.;
  nodes(1) = {1., 0., 0.};
  nodes(2) = {1., 1., 0.};
  nodes(3) = {1., 0., 0.};
  hexed::Array<double> tangent_averages({4, 3});
  tangent_averages(0) = .1;
  tangent_averages(1) = .2;
  tangent_averages(2) = .3;
  tangent_averages(3) = .2;
  hexed::Array<double> tangent_radii = hexed::Array<double>::make(.4, .5, .6, .5);
  std::vector<std::shared_ptr<hexed::Geom_edge>> edges;
  edges.push_back(std::make_shared<hexed::Tree_curve_edge>(nodes(0, 2), tangent_averages(0, 2), tangent_radii(0, 2)));
  edges.push_back(std::make_shared<hexed::Tree_curve_edge>(nodes(2, 4), tangent_averages(2, 4), tangent_radii(2, 4)));
  hexed::Compound_edge edge(edges, {false, true});
  REQUIRE(edge.arc_length(0.) == Catch::Approx(0.).scale(1.));
  REQUIRE(edge.arc_length(.25) == Catch::Approx(.5).scale(1.));
  REQUIRE(edge.arc_length(.75) == Catch::Approx(1.5).scale(1.));
  REQUIRE(edge.arc_length(1.) == Catch::Approx(2.).scale(1.));
  REQUIRE((edge.point(0.) - hexed::Mat<3>::Zero()).norm() < 1e-6);
  REQUIRE((edge.point(.25) - hexed::Mat<3>{.5, 0., 0.}).norm() < 1e-6);
  REQUIRE((edge.point(.75) - hexed::Mat<3>{1., .5, 0.}).norm() < 1e-6);
  REQUIRE((edge.point(1.) - hexed::Mat<3>{1., 1., 0.}).norm() < 1e-6);
  REQUIRE((edge.tangent_average(.25) - hexed::Mat<3>{.15, .15, .15}).norm() < 1e-6);
  REQUIRE(edge.tangent_radius(1.) == Catch::Approx(.6).scale(1.));
  REQUIRE(edge.arg_nearest_point(hexed::Mat<3>{0.2, 0.01, 0.}) == Catch::Approx(.1));
  REQUIRE(edge.arg_nearest_point(hexed::Mat<3>{1.01, .2, 0.}) == Catch::Approx(.6));
}

TEST_CASE("Compound_geom") {
  std::vector<hexed::Surface_geom*> geoms;
  geoms.push_back(new hexed::Hypersphere(Eigen::Vector3d::Zero(), 1.));
  geoms.push_back(new hexed::Hypersphere(Eigen::Vector3d::Unit(0), 1.));
  hexed::Compound_geom geom(geoms);
  REQUIRE_THAT(geom.nearest_point(Eigen::Vector3d{-2., 0., 0.}).point(),
               Catch::Matchers::RangeEquals(Eigen::Vector3d{-1., 0., 0.}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE_THAT(geom.nearest_point(Eigen::Vector3d{3., 0., 0.}).point(),
               Catch::Matchers::RangeEquals(Eigen::Vector3d{2., 0., 0.}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE(geom.nearest_point(Eigen::Vector3d{3., 0., 0.}, 0.1).empty());
  REQUIRE_THAT(geom.intersections(Eigen::Vector3d::Zero(), Eigen::Vector3d::Unit(0)),
               Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-1., 0., 1., 2.},
               hexed::math::Approx_equal(0., 1e-12)));
}

TEST_CASE("Hypersphere") {
  hexed::Hypersphere hype(Eigen::Vector2d{.1, .2}, .5);
  REQUIRE_THAT(hype.nearest_point(Eigen::Vector2d{-.5, -.6}).point(),
               Catch::Matchers::RangeEquals(Eigen::Vector2d{-.2, -.2}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE(hype.nearest_point(Eigen::Vector2d{-.5, -.6}, 0.1).empty());
  // intersections: {.4, .6}, {.5, .5}. diff: {.1, -.1}
  REQUIRE_THAT(hype.intersections(Eigen::Vector2d{.6, .4}, Eigen::Vector2d{.7, .3}),
               Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-1., -2.},
               hexed::math::Approx_equal(0., 1e-12)));
}
