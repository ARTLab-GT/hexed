#include <catch2/catch_all.hpp>
#include <hexed/Geom_edge.hpp>

TEST_CASE("Geom_edge") {
  REQUIRE_THROWS(hexed::Geom_edge(hexed::Array<double>({4, 5})));
  REQUIRE_THROWS(hexed::Geom_edge(hexed::Array<double>({2, 2, 3})));
  REQUIRE_THROWS(hexed::Geom_edge(hexed::Array<double>({1, 3})));
  std::vector<double> points {1., 1., 1., 2., 2., 2., 3., 3., 2.};

  SECTION("2 points") {
    hexed::Geom_edge edge(hexed::Array<double>({2, 3}, points.data()));
    std::vector<double> correct {1., 1., 1., 2., 2., 2.};
    REQUIRE_THAT(edge.points(), Catch::Matchers::RangeEquals(correct, hexed::math::Approx_equal()));
    correct = std::vector<double>{0., std::sqrt(3.)};
    REQUIRE_THAT(edge.arc_len(), Catch::Matchers::RangeEquals(correct, hexed::math::Approx_equal(0., 1e-8)));
  }

  SECTION("3 points") {
    hexed::Geom_edge edge(hexed::Array<double>({3, 3}, points.data()));
    REQUIRE_THAT(edge.points(), Catch::Matchers::RangeEquals(points, hexed::math::Approx_equal()));
    std::vector<double> correct {0., std::sqrt(3.), std::sqrt(3.) + std::sqrt(2.)};
    REQUIRE_THAT(edge.arc_len(), Catch::Matchers::RangeEquals(correct, hexed::math::Approx_equal(0., 1e-8)));
    hexed::Mat<3> to{2.1, 2.1, 2.1};
    auto node = edge.nearest(to);
    REQUIRE_THAT(node.pos, Catch::Matchers::RangeEquals(std::vector<double>{2., 2., 2.}, hexed::math::Approx_equal()));
    REQUIRE(node.arc_len == Catch::Approx(std::sqrt(3.)));
    REQUIRE(node.index == 1);
    node = edge.nearest(to, std::sqrt(3.) + .01);
    REQUIRE(node.arc_len == Catch::Approx(std::sqrt(3.) + std::sqrt(2.)));
    REQUIRE(node.index == 2);
    node = edge.nearest(to, 0., .01);
    REQUIRE(node.arc_len == Catch::Approx(0.).margin(1e-6));
    REQUIRE(node.index == 0);
  }
}
