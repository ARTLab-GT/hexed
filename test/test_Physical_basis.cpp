#include <catch2/catch.hpp>
#include <Physical_basis.hpp>

TEST_CASE("Physical_basis")
{
  SECTION("size checking")
  {
    std::vector<double> points (50);
    std::vector<cartdg::Physical_basis> bases;
    REQUIRE_THROWS(bases.emplace_back(1, 5, points));
    REQUIRE_THROWS(bases.emplace_back(2, 4, points));
    REQUIRE_THROWS(bases.emplace_back(2, 6, points));
    bases.emplace_back(2, 5, points);
    std::vector<double> short_points (4*16);
    REQUIRE_THROWS(bases.emplace_back(4, 2, short_points));
  }
}
