#include <catch2/catch_all.hpp>
#include <hexed/Iges_parser.hpp>

TEST_CASE("Iges_parser") {
  for (std::string file_name : {"default_delims", "custom_param", "custom_record", "custom_both"}) {
    hexed::Iges_parser file("../test_assets/" + file_name + ".iges");
    auto start_sec = file.section(hexed::Iges_parser::start);
    REQUIRE(start_sec.size() == 2);
    REQUIRE(start_sec[0].size() == 1);
    REQUIRE(start_sec[1].size() == 1);
    REQUIRE(start_sec[0][0] == "first start entry                                                       ");
    REQUIRE(start_sec[1][0] == "second start entry                                                      ");
    REQUIRE(&file.entry(hexed::Iges_parser::start, 1) == &start_sec[0]);
    REQUIRE(&file.entry(hexed::Iges_parser::start, 2) == &start_sec[1]);
    REQUIRE(file.section(hexed::Iges_parser::global).size() == 1);
    REQUIRE_THAT(file.entry(hexed::Iges_parser::global, 1), Catch::Matchers::RangeEquals(std::vector<std::string> {
      "12Hfirst,_entry",
      "12Hsecond entry",
      "5Hthird",
      "25",
      "",
      "9.8",
    }));
    if (file_name == "default_delims") {
      auto param_sec = file.section(hexed::Iges_parser::parameter);
      REQUIRE(param_sec.size() == 2);
      REQUIRE(&file.entry(hexed::Iges_parser::parameter, 1) == &param_sec[0]);
      REQUIRE(&file.entry(hexed::Iges_parser::parameter, 3) == &param_sec[1]);
      REQUIRE_THAT(param_sec[0], Catch::Matchers::RangeEquals(std::vector<std::string> {"0", "1.2"}));
      REQUIRE_THAT(param_sec[1], Catch::Matchers::RangeEquals(std::vector<std::string> {"0"}));
    }
  }
}
