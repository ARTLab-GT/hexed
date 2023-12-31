#include <catch2/catch_all.hpp>
#include <hexed/Interpreter.hpp>

TEST_CASE("Interpreter")
{
  hexed::Interpreter inter(std::vector<std::string>{});
  inter.exec("\n");
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
  inter.exec("str = {\\{}");
  REQUIRE(inter.variables->lookup<std::string>("str").value() == "{");
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
  inter.exec("code = -1.1 + {name} + 3");
  REQUIRE(inter.variables->lookup<std::string>("code").value() == "-1.1name3");
  inter.exec("result = 2*3 + 1*2 - -3*-3*3");
  REQUIRE(inter.variables->lookup<int>("result").value() == -19);
  inter.exec("result = (1 + 2*2)*2");
  REQUIRE(inter.variables->lookup<int>("result").value() == 10);
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
  inter.exec("result = 0; exit; result = 1");
  REQUIRE(inter.variables->lookup<int>("result").value() == 0);
  inter.exec("size = #{prandtl}");
  REQUIRE(inter.variables->lookup<int>("size").value() == 7);
  inter.exec("char = {prandtl}#1");
  REQUIRE(inter.variables->lookup<std::string>("char").value() == "r");
  // test that `exec` calls are properly mutex'd
  #pragma omp parallel for
  for (int i = 0; i < 10; ++i) inter.exec("result = " + std::to_string(i));
  int result = inter.variables->lookup<int>("result").value();
  REQUIRE(result >= 0);
  REQUIRE(result < 10);
  inter.exec("except = {result = 21}");
  inter.exec("result = nonexistant");
  REQUIRE(inter.variables->lookup<int>("result") == 21);
  // test evaluation of assignments
  inter.exec("a = b = 2");
  REQUIRE(inter.variables->lookup<int>("a") == 2);
  inter.exec("a = (b = 2) + 2");
  REQUIRE(inter.variables->lookup<int>("a") == 4);
  inter.exec("a = (b = 3; c = 5)");
  REQUIRE(inter.variables->lookup<int>("a") == 5);
  inter.exec("a"); // expressions don't have to contain assignments
  inter.exec("x = (0; {})");

  SECTION("standard library")
  {
    hexed::Interpreter test;
    test.exec("$(read {../test/test_builtin.hil})");
    REQUIRE(test.variables->lookup<int>("triangle").value() == 10);
    REQUIRE(test.variables->lookup<int>("cube").value() == 27);
  }
}
