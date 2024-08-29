#include <catch2/catch_all.hpp>
#include <hexed/Tree_curve.hpp>

TEST_CASE("Tree_curve") {
  REQUIRE_THROWS(hexed::Tree_curve(hexed::Array<double>({66, 3})));
  hexed::Array<double> nodes({65, 3});
  for (int i = 0; i < 65; ++i) nodes(i).vector() << i/64., i*i/64./64., 3.;
  hexed::Tree_curve curve(nodes, 2);
  REQUIRE(curve.skip_levels() == 2);
  REQUIRE(curve.nodes().shape()[0] == 65);
  REQUIRE(curve.nodes().shape()[1] == 3);
  REQUIRE(curve.segments(3)[0].nodes.shape()[0] == 9);
  REQUIRE(curve.segments(3)[0].nodes.shape()[1] == 3);
  REQUIRE(curve.segments(3)[7].nodes.shape()[0] == 9);
  REQUIRE(curve.segments(3)[0].nodes(0).data() == curve.root().nodes(0).data());
  REQUIRE(curve.segments(3)[7].nodes(8).data() == curve.root().nodes(64).data());
  REQUIRE(curve.root().segments[0].segments[1].nodes(0).data() == nodes(16).data());
  for (int level = 0; level < 4; ++level) {
    for (int i_segment = 0; i_segment < hexed::math::pow(2, level); ++i_segment) {
      auto& seg = curve.segments(level)[i_segment];
      for (int i_node = 0; i_node < hexed::math::pow(2, 6 - level) + 1; ++i_node) {
        REQUIRE((seg.nodes(i_node).vector() - seg.center).norm() < seg.radius + 1e-8);
      }
    }
  }
  REQUIRE(curve.segments(3)[0].segments.size() == 0);

  hexed::Array<double> circle_nodes({1025, 3});
  for (int i = 0; i < 1025; ++i) {
    double angle = i*2.*M_PI/1025;
    circle_nodes(i)[0] = std::cos(angle);
    circle_nodes(i)[1] = std::sin(angle);
    circle_nodes(i)[2] = 1.;
  }
  hexed::Tree_curve circle(circle_nodes);
  hexed::Mat<3> point {std::sqrt(.5) + .1, std::sqrt(.5) + .1, 1.1};
  REQUIRE_THAT(circle.nearest_point(point, .3).point(),
               Catch::Matchers::RangeEquals(hexed::Mat<3>{std::sqrt(.5), std::sqrt(.5), 1.},
                                            hexed::math::Approx_equal(0., 2e-3)));
  REQUIRE(circle.nearest_point(point, .1).empty());
}
