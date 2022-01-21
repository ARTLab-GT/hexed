#include <catch2/catch.hpp>
#include <Physical_basis.hpp>

TEST_CASE("Physical_basis")
{
  SECTION("input size checking")
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

  SECTION("basis size calculation")
  {
    {
      std::vector<double> points(2);
      cartdg::Physical_basis basis {2, 1, points};
      REQUIRE(basis.size() == 1);
    }
    {
      std::vector<double> points(8);
      cartdg::Physical_basis basis {2, 2, points};
      REQUIRE(basis.size() == 3);
    }
    {
      std::vector<double> points(32);
      cartdg::Physical_basis basis {2, 4, points};
      REQUIRE(basis.size() == 10);
    }
    {
      std::vector<double> points(3);
      cartdg::Physical_basis basis {3, 1, points};
      REQUIRE(basis.size() == 1);
    }
    {
      std::vector<double> points(24);
      cartdg::Physical_basis basis {3, 2, points};
      REQUIRE(basis.size() == 4);
    }
    {
      std::vector<double> points(192);
      cartdg::Physical_basis basis {3, 4, points};
      REQUIRE(basis.size() == 20);
    }
  }
}
