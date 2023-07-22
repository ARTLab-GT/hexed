#include <catch2/catch_all.hpp>
#include <hexed/Surface_geom.hpp>

TEST_CASE("Compound_geom")
{
  std::vector<hexed::Surface_geom*> geoms;
  geoms.push_back(new hexed::Hypersphere(Eigen::Vector3d::Zero(), 1.));
  geoms.push_back(new hexed::Hypersphere(Eigen::Vector3d::Unit(0), 1.));
  hexed::Compound_geom geom(geoms);
  REQUIRE_THAT(geom.nearest_point(Eigen::Vector3d{-2., 0., 0.}),
               Catch::Matchers::RangeEquals(Eigen::Vector3d{-1., 0., 0.}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE_THAT(geom.nearest_point(Eigen::Vector3d{3., 0., 0.}),
               Catch::Matchers::RangeEquals(Eigen::Vector3d{2., 0., 0.}, hexed::math::Approx_equal(0., 1e-12)));
  REQUIRE_THAT(geom.intersections(Eigen::Vector3d::Zero(), Eigen::Vector3d::Unit(0)),
               Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-1., 0., 1., 2.}, hexed::math::Approx_equal(0., 1e-12)));
}

TEST_CASE("Hypersphere")
{
  hexed::Hypersphere hype(Eigen::Vector2d{.1, .2}, .5);
  REQUIRE_THAT(hype.nearest_point(Eigen::Vector2d{-.5, -.6}),
               Catch::Matchers::RangeEquals(Eigen::Vector2d{-.2, -.2}, hexed::math::Approx_equal(0., 1e-12)));
  // intersections: {.4, .6}, {.5, .5}. diff: {.1, -.1}
  REQUIRE_THAT(hype.intersections(Eigen::Vector2d{.6, .4}, Eigen::Vector2d{.7, .3}),
               Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-1., -2.}, hexed::math::Approx_equal(0., 1e-12)));
}
