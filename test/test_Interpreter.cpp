#include <catch2/catch_all.hpp>
#include <hexed/Interpreter.hpp>
#include <iostream>


TEST_CASE("Interpreter")
{
  hexed::Interpreter inter;
  inter.exec("shock_wave = 7\n  boundary0layer=14;interaction = 1.2");
  REQUIRE(inter.variables->lookup<int>("shock_wave").value() == 7);
  REQUIRE(inter.variables->lookup<int>("boundary0layer").value() == 14);
  REQUIRE(inter.variables->lookup<double>("interaction").value() == Catch::Approx(1.2));
  inter.exec("interaction = .2");
  REQUIRE(inter.variables->lookup<double>("interaction").value() == Catch::Approx(.2));
  inter.exec("interaction = 1e-3");
  REQUIRE(inter.variables->lookup<double>("interaction").value() == Catch::Approx(1e-3));
  inter.exec("ludwig = {prandtl}\ntitle = {consider\n{Phlebas}!\\}\\\\}");
  REQUIRE(inter.variables->lookup<std::string>("ludwig").value() == "prandtl");
  REQUIRE(inter.variables->lookup<std::string>("title").value() == "consider\n{Phlebas}!}\\");
  inter.exec("boundary0layer = shock_wave");
  REQUIRE(inter.variables->lookup<int>("boundary0layer").value() == 7);
  inter.exec("boundary0layer = -shock_wave");
  REQUIRE(inter.variables->lookup<int>("boundary0layer").value() == -7);
  inter.exec("interaction = 0.1*boundary0layer");
  REQUIRE(inter.variables->lookup<double>("interaction").value() == Catch::Approx(-.7));
  inter.exec("result = 2*7*-9/3");
  REQUIRE(inter.variables->lookup<int>("result").value() == -42);
  inter.exec("result = 1 + 1");
  REQUIRE(inter.variables->lookup<int>("result").value() == 2);
  inter.exec("result = 12%7");
  REQUIRE(inter.variables->lookup<int>("result").value() == 5);
  inter.exec("ludwig = ludwig + {-glauert}");
  REQUIRE(inter.variables->lookup<std::string>("ludwig").value() == "prandtl-glauert");
  inter.exec("code = -1.0 + {name} + 3");
  REQUIRE(inter.variables->lookup<std::string>("code").value() == "-1.000000name3");
  inter.exec("result = 2*3 + 1*2 - -3*-3*3");
  REQUIRE(inter.variables->lookup<int>("result").value() == -19);
  inter.exec("result = (1 + 2)*2");
  REQUIRE(inter.variables->lookup<int>("result").value() == 6);
  inter.exec("result = 6/3*2");
  REQUIRE(inter.variables->lookup<int>("result").value() == 4);
  inter.exec("result = !0");
  REQUIRE(inter.variables->lookup<int>("result").value() == 1);
  inter.exec("result = !2");
  REQUIRE(inter.variables->lookup<int>("result").value() == 0);
  inter.exec("root = sqrt 9");
  REQUIRE(inter.variables->lookup<double>("root").value() == Catch::Approx(3.));
  inter.exec("root = sqrt 4.");
  REQUIRE(inter.variables->lookup<double>("root").value() == Catch::Approx(2.));
  inter.exec("code = {2 + 7}");
  inter.exec("result = $code");
  REQUIRE(inter.variables->lookup<int>("result").value() == 9);
  inter.exec("code = {result = result + 1}");
  inter.exec("$code");
  REQUIRE(inter.variables->lookup<int>("result").value() == 10);
  inter.exec("result = 1 + ${${1 + 1}}");
  REQUIRE(inter.variables->lookup<int>("result").value() == 3);
  inter.exec("= 1 + 1  ");
  inter.exec("result = 6; result = result*result");
  REQUIRE(inter.variables->lookup<int>("result").value() == 36);
  inter.exec("result = -7.1 < 5");
  REQUIRE(inter.variables->lookup<int>("result").value() == 1);
  inter.exec("result = 3 == 5");
  REQUIRE(inter.variables->lookup<int>("result").value() == 0);
  inter.exec("result = 3 == 2 + 1");
  REQUIRE(inter.variables->lookup<int>("result").value() == 1);
  inter.exec("result = 3 >= 5");
  REQUIRE(inter.variables->lookup<int>("result").value() == 0);
  inter.exec("result = 1 & 0");
  REQUIRE(inter.variables->lookup<int>("result").value() == 0);
  inter.exec("result = 1 | 0");
  REQUIRE(inter.variables->lookup<int>("result").value() == 1);
  inter.exec("result = 2^5");
  REQUIRE(inter.variables->lookup<int>("result").value() == 32);
  inter.exec("real = 1e-4^0.25");
  REQUIRE(inter.variables->lookup<double>("real").value() == Catch::Approx(0.1));
  inter.exec("result = {prandtl} == {Prandtl}");
  REQUIRE(inter.variables->lookup<int>("result").value() == 0);
  inter.exec("result = {prandtl} == {prandtl}");
  REQUIRE(inter.variables->lookup<int>("result").value() == 1);
  REQUIRE_THROWS(inter.exec("$read {non_existant.hil}"));
  inter.exec("result = 0; =exit; result = 1");
  REQUIRE(inter.variables->lookup<int>("result").value() == 0);

  SECTION("standard library")
  {
    hexed::Interpreter test;
    test.exec("$(read {../test/test_std.hil})");
    REQUIRE(test.variables->lookup<std::string>("result0").value() == "yes");
    REQUIRE(test.variables->lookup<std::string>("result1").value() == "yesyes");
    REQUIRE(test.variables->lookup<int>("triangle").value() == 10);
    REQUIRE(test.variables->lookup<int>("cube").value() == 27);
    test.exec("$repl");
  }
}
