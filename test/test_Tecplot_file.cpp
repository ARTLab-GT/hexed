#include <catch2/catch.hpp>
#include <Tecplot_file.hpp>

TEST_CASE("Tecplot_file")
{
  {
    cartdg::Tecplot_file file0 {"test0", 3, 5, 4, 0.2};
    REQUIRE_THROWS(cartdg::Tecplot_file("test0", 3, 5, 4, 0.2));
  }
  cartdg::Tecplot_file file0 {"test0", 3, 5, 4, 0.2};
}
