#include <catch2/catch.hpp>

#include <static_math.hpp>

TEST_CASE("static_pow")
{
  {
    constexpr int result = cartdg::static_math::pow(2, 0);
    REQUIRE(result == 1);
  }
  {
    constexpr int result = cartdg::static_math::pow(2, 1);
    REQUIRE(result == 2);
  }
  {
    constexpr int result = cartdg::static_math::pow(2, -1);
    REQUIRE(result == 0);
  }
  {
    constexpr int result = cartdg::static_math::pow(3, 4);
    REQUIRE(result == 81);
  }
  {
    constexpr double result = cartdg::static_math::pow(1.5, 2);
    REQUIRE(result == 2.25);
  }
  {
    constexpr double result = cartdg::static_math::pow(2., -1);
    REQUIRE(result == 0.5);
  }
}
