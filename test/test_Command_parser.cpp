#include <catch2/catch_all.hpp>
#include <hexed/Command_parser.hpp>
#include <iostream>


TEST_CASE("Command_parser")
{
  hexed::Command_parser parser;
  parser.exec("shock_wave = 7");
  REQUIRE(parser.variables->lookup<int>("shock_wave").value() == 7);
}
