#include <catch2/catch_all.hpp>
#include <hexed/Array.hpp>

TEST_CASE("Array")
{
  hexed::Array<double, true> arr0({2, 3, 4});
  REQUIRE(arr0.order() == 3);
  REQUIRE_THAT(arr0.shape(), Catch::Matchers::RangeEquals(std::vector<int>{2, 3, 4}));
  REQUIRE(arr0.size() == 24);

  hexed::Array<double> zero_size({});
  REQUIRE(zero_size.shape().empty());
  REQUIRE(zero_size.size() == 0);
  REQUIRE(arr0.data()[0] == 0.);
  REQUIRE(arr0[23] == 0.);

  auto view = arr0();
  view[1] = 0.6;
  REQUIRE(arr0[1] == Catch::Approx(0.6));
  arr0[20] = -1.;
  REQUIRE(arr0(0).order() == 2);
  REQUIRE(arr0(1)(0).size() == 4);
  REQUIRE_THAT(arr0(1).shape(), Catch::Matchers::RangeEquals(std::vector<int>{3, 4}));
  REQUIRE(arr0(1)(2)[0] == Catch::Approx(-1.));
  REQUIRE_THROWS(arr0(0)(0)(0)(0));
  REQUIRE_THROWS(arr0(2));
  REQUIRE_THROWS(arr0(0)(0)(0)[0]);
  REQUIRE_THROWS(arr0[24]);
}
