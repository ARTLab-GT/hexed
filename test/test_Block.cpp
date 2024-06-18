#include <catch2/catch_all.hpp>
#include <hexed/Block.hpp>
#include <hexed/config.hpp>
#include <hexed/Equidistant.hpp>

TEST_CASE("Block")
{
  static_assert(hexed::config::max_row_size >= 4); // these tests require a row size of at least 4
  hexed::next::Vertex vert0({0.1, -.3, .2}, 4);
  REQUIRE(vert0.n_dim == 0);
  REQUIRE(vert0.row_size == 4);
  REQUIRE_THAT(vert0.point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, -.3, .2}, hexed::math::Approx_equal()));
  hexed::next::Vertex vert1({.3, -.1, .4}, 4);
  hexed::next::Edge edge0(vert0, vert1, std::make_shared<hexed::Equidistant>(4));
  REQUIRE(edge0.n_dim == 1);
  REQUIRE(edge0.row_size == 4);
  auto test_interp = [&](){
    for (int i = 0; i < 4; ++i) {
      REQUIRE_THAT(edge0.point({i}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, -.3, .2} + i*hexed::Mat<3>::Constant(.2/3.), hexed::math::Approx_equal()));
    }
  };
  test_interp();
}
