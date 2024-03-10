#include <catch2/catch_all.hpp>
#include <hexed/Array.hpp>

TEST_CASE("Array")
{
  hexed::Array<double> arr0({2, 3, 4});
  REQUIRE(arr0.order() == 3);
  REQUIRE_THAT(arr0.shape(), Catch::Matchers::RangeEquals(std::vector<int>{2, 3, 4}));
  REQUIRE(arr0.size() == 24);
  hexed::Array<double> zero_size({});
  REQUIRE(zero_size.shape().empty());
  REQUIRE(zero_size.size() == 0);
}
