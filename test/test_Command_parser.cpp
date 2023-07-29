#include <catch2/catch_all.hpp>
#include <hexed/Command_parser.hpp>
#include <iostream>


TEST_CASE("Command_parser")
{
  hexed::Command_parser parser;
  parser.exec("shock_wave = 7\n  boundary0layer=14 \n\ninteraction = 1.2");
  REQUIRE(parser.variables->lookup<int>("shock_wave").value() == 7);
  REQUIRE(parser.variables->lookup<int>("boundary0layer").value() == 14);
  REQUIRE(parser.variables->lookup<double>("interaction").value() == Catch::Approx(1.2));
  parser.exec("interaction = .2");
  REQUIRE(parser.variables->lookup<double>("interaction").value() == Catch::Approx(.2));
  parser.exec("interaction = 1e3");
  REQUIRE(parser.variables->lookup<double>("interaction").value() == Catch::Approx(1e3));
  parser.exec("ludwig = \"prandtl\"\ntitle = \"consider\nPlebas\"\"!\"\"\"");
  REQUIRE(parser.variables->lookup<std::string>("ludwig").value() == "prandtl");
  REQUIRE(parser.variables->lookup<std::string>("title").value() == "consider\nPlebas\"!\"");
  parser.exec("boundary0layer = shock_wave");
  REQUIRE(parser.variables->lookup<int>("boundary0layer").value() == 7);
}
