#include <catch2/catch_all.hpp>
#include <hexed/utils.hpp>

TEST_CASE("format_str") {
  REQUIRE(hexed::format_str(100, "%.2f == %d", M_PI, 3) == std::string{"3.14 == 3"});
  REQUIRE_THROWS(hexed::format_str(3, "g == %d", 10)); // too many characters
  REQUIRE(hexed::file_extension("archive.tar.gz") == "gz");
  REQUIRE(hexed::file_extension("model.STL") == "stl");
}
