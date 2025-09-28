#include <catch2/catch_all.hpp>
#include <hexed/utils.hpp>
#include <hexed/math.hpp>

TEST_CASE("file_extension") {
  REQUIRE(hexed::file_extension("archive.tar.gz") == "gz");
  REQUIRE(hexed::file_extension("model.STL") == "stl");
}

TEST_CASE("to_lower") {
  REQUIRE(hexed::to_lower("SOmE cHarac-teRs.") == "some charac-ters.");
  REQUIRE(hexed::to_lower("") == "");
  REQUIRE(hexed::to_lower("\n") == "\n");
}

TEST_CASE("to_mat") {
  std::vector<double> vec {.1, -.3, .2};
  REQUIRE_THAT(hexed::to_mat(vec), Catch::Matchers::RangeEquals(vec, hexed::math::Approx_equal()));
}

TEST_CASE("resize") {
  hexed::Mat<4> vec {.2, -.1, .03, 6.};
  REQUIRE_THAT(hexed::resize(vec, 2),
               Catch::Matchers::RangeEquals(std::vector<double>{.2, -.1}, hexed::math::Approx_equal(0., 1e-10)));
  REQUIRE_THAT(hexed::resize(vec, 6),
               Catch::Matchers::RangeEquals(std::vector<double>{.2, -.1, .03, 6., 0., 0.},
                                            hexed::math::Approx_equal(0., 1e-10)));
}

TEST_CASE("str_cat") {
  REQUIRE(hexed::str_cat(42) == "42");
  REQUIRE(hexed::str_cat(8314.0, "ludwig prandtl", true) == "+8.314000e+03ludwig prandtltrue");
}
