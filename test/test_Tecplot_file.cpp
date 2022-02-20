#include <catch2/catch.hpp>
#include <Tecplot_file.hpp>

TEST_CASE("Tecplot_file")
{
  double pos = 0.;
  double state = 0.;
  std::vector<std::string> var_names {"var0", "var1"};
  {
    cartdg::Tecplot_file file0 {"test0", 3, var_names, 0.2};
    REQUIRE_THROWS(cartdg::Tecplot_file("test0", 3, var_names, 0.2)); // only one `Tecplot_file` can exist at a time.
    {
      cartdg::Tecplot_file::Structured_block zone (file0, 1);
      REQUIRE_THROWS(cartdg::Tecplot_file::Structured_block(file0, 1));
      zone.write(&pos, &state);
    }
  }
  cartdg::Tecplot_file file0 {"test0", 3, var_names, 0.2};
  cartdg::Tecplot_file::Structured_block zone (file0, 1);
  zone.write(&pos, &state);
}
