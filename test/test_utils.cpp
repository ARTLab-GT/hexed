#include <catch2/catch.hpp>
#include <hexed/utils.hpp>

TEST_CASE("format_str")
{
  REQUIRE(hexed::format_str(100, "%.2f == %d", M_PI, 3) == std::string{"3.14 == 3"});
  REQUIRE_THROWS(hexed::format_str(3, "g == %d", 10)); // too many characters
}
