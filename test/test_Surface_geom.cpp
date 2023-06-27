#include <catch2/catch_all.hpp>
#include <hexed/Surface_geom.hpp>

TEST_CASE("Hypersphere")
{
  hexed::Hypersphere hype(Eigen::Vector2d{.1, .2}, .5);
  REQUIRE_THAT(hype.nearest_point(Eigen::Vector2d{-.5, -.6}),
               Catch::Matchers::RangeEquals(Eigen::Vector2d{-.2, -.2}, hexed::math::Approx_equal(1e-12)));
  // intersections: {.4, .6}, {.5, .5}. diff: {.1, -.1}
  REQUIRE_THAT(hype.intersections(Eigen::Vector2d{.6, .4}, Eigen::Vector2d{.7, .3}),
               Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-1., -2.}, hexed::math::Approx_equal(1e-12)));
}
