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
  parser.exec("interaction = 1e-3");
  REQUIRE(parser.variables->lookup<double>("interaction").value() == Catch::Approx(1e-3));
  parser.exec("ludwig = \"prandtl\"\ntitle = \"consider\nPhlebas\"\"!\"\"\"");
  REQUIRE(parser.variables->lookup<std::string>("ludwig").value() == "prandtl");
  REQUIRE(parser.variables->lookup<std::string>("title").value() == "consider\nPhlebas\"!\"");
  parser.exec("boundary0layer = shock_wave");
  REQUIRE(parser.variables->lookup<int>("boundary0layer").value() == 7);
  parser.exec("boundary0layer = -shock_wave");
  REQUIRE(parser.variables->lookup<int>("boundary0layer").value() == -7);
  parser.exec("interaction = 0.1*boundary0layer");
  REQUIRE(parser.variables->lookup<double>("interaction").value() == Catch::Approx(-.7));
  parser.exec("result = 2*7*-9/3");
  REQUIRE(parser.variables->lookup<int>("result").value() == -42);
  parser.exec("result = 1 + 1");
  REQUIRE(parser.variables->lookup<int>("result").value() == 2);
  parser.exec("result = 2*3 + 1*2");
  REQUIRE(parser.variables->lookup<int>("result").value() == 8);
}
