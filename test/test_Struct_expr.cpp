#include <catch2/catch_all.hpp>
#include <hexed/Struct_expr.hpp>

TEST_CASE("Struct_expr")
{
  hexed::Interpreter inter;
  inter.exec("macro = {var}");
  REQUIRE_THROWS(hexed::Struct_expr(inter, "$macro = 0"));
  hexed::Struct_expr expr(inter, "var0 = 2; var1 = 2 + var0\n\nvar2 =\t${10}");
  REQUIRE_THAT(expr.names, Catch::Matchers::RangeEquals(std::vector<std::string>{"var0", "var1", "var2"}));
  REQUIRE_THAT(expr.exprs, Catch::Matchers::RangeEquals(std::vector<std::string>{"2", "2 + var0", "${10}"}));
  REQUIRE_THAT(expr.eval(), Catch::Matchers::RangeEquals(std::vector<double>{2, 4, 10}, hexed::math::Approx_equal()));
}
