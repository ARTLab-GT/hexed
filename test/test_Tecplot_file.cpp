#include <hexed/config.hpp>
#if HEXED_USE_TECPLOT
#include <catch2/catch_all.hpp>
#include <hexed/Tecplot_file.hpp>

TEST_CASE("Tecplot_file")
{
  double pos = 0.;
  double state = 0.;
  std::vector<std::string> var_names {"var0", "var1"};
  {
    hexed::Tecplot_file file0 {"test0", 3, var_names, 0.2};
    hexed::Tecplot_file::Structured_block zone (file0, 1);
    zone.write(&pos, &state);
  }
  hexed::Tecplot_file file0 {"test0", 3, var_names, 0.2};
  hexed::Tecplot_file::Structured_block zone (file0, 1);
  zone.write(&pos, &state);
}
#endif
