#include <catch2/catch_all.hpp>
#include <hexed/Interpreter.hpp>
#include <iostream>


TEST_CASE("Interpreter")
{
  hexed::Interpreter inter;
  inter.exec("shock_wave = 7\n  boundary0layer=14 \n\ninteraction = 1.2");
  REQUIRE(inter.variables->lookup<int>("shock_wave").value() == 7);
  REQUIRE(inter.variables->lookup<int>("boundary0layer").value() == 14);
  REQUIRE(inter.variables->lookup<double>("interaction").value() == Catch::Approx(1.2));
  inter.exec("interaction = .2");
  REQUIRE(inter.variables->lookup<double>("interaction").value() == Catch::Approx(.2));
  inter.exec("interaction = 1e-3");
  REQUIRE(inter.variables->lookup<double>("interaction").value() == Catch::Approx(1e-3));
  inter.exec("ludwig = \"prandtl\"\ntitle = \"consider\nPhlebas\"\"!\"\"\"");
  REQUIRE(inter.variables->lookup<std::string>("ludwig").value() == "prandtl");
  REQUIRE(inter.variables->lookup<std::string>("title").value() == "consider\nPhlebas\"!\"");
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
  inter.exec("code = \"2 + 7\"");
  inter.exec("result = $code");
  REQUIRE(inter.variables->lookup<int>("result").value() == 9);
  inter.exec("code = \"result = result + 1\"");
  inter.exec("$code");
  REQUIRE(inter.variables->lookup<int>("result").value() == 10);

  inter.exec("$(read \"../include/std.hil\")");
  REQUIRE(inter.variables->lookup<int>("sum").value() == 2);
}
