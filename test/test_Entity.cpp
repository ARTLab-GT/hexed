#include <catch2/catch_all.hpp>
#include <hexed/Entity.hpp>
#include <hexed/Gauss_lobatto.hpp>

TEST_CASE("Entity")
{
  // this test requires a max row size of at least 3
  static_assert(hexed::config::max_row_size >= 3);
  hexed::Cartesian ent3(3, 3, std::make_shared<hexed::Gauss_lobatto>(3), hexed::Array<int>::make(1, -3, -1), .1);
  REQUIRE_THAT(ent3.nominal_position().vector(), Catch::Matchers::RangeEquals(std::vector<int>{1, -3, -1}));
  auto points = ent3.points();
  REQUIRE_THAT(points.shape(), Catch::Matchers::RangeEquals(std::vector<int>{3, 3, 3, 3}));

  REQUIRE(points(0)(0)(0)[0] == Catch::Approx(.1));
  REQUIRE(points(1)(0)(0)[0] == Catch::Approx(-.3));
  REQUIRE(points(2)(0)(0)[0] == Catch::Approx(-.1));

  REQUIRE(points(0)(2)(0)[0] == Catch::Approx(.2));
  REQUIRE(points(1)(2)(0)[0] == Catch::Approx(-.3));
  REQUIRE(points(2)(2)(0)[0] == Catch::Approx(-.1));

  REQUIRE(points(0)(0)(2)[0] == Catch::Approx(.1));
  REQUIRE(points(1)(0)(2)[0] == Catch::Approx(-.2));
  REQUIRE(points(2)(0)(2)[0] == Catch::Approx(-.1));

  REQUIRE(points(0)(1)(1)[1] == Catch::Approx(.15));
  REQUIRE(points(1)(1)(1)[1] == Catch::Approx(-.25));
  REQUIRE(points(2)(1)(1)[1] == Catch::Approx(-.05));
}
