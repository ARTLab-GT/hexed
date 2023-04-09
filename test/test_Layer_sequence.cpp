#include <catch2/catch_all.hpp>
#include <hexed/Layer_sequence.hpp>

TEST_CASE("Layer_sequence")
{
  hexed::Layer_sequence ls(1e-6, 4);
  // test total height is 1
  REQUIRE(ls.cumulative_height(ls.n_layers()) == Catch::Approx(1.).epsilon(1e-10));
  // test near-wall layers
  for (int i = 0; i < 4; ++i) REQUIRE(ls.spacing(i) == Catch::Approx(1e-6));
  REQUIRE(ls.cumulative_height(0) == Catch::Approx(0.).scale(1e-6));
  REQUIRE(ls.cumulative_height(4) == Catch::Approx(4e-6));
  // test growth ratios
  REQUIRE(ls.spacing(ls.n_layers() - 1) > .49);
  for (int i = 1; i < ls.n_layers(); ++i) {
    REQUIRE(ls.spacing(i)/ls.spacing(i - 1) < 2.01);
    REQUIRE(ls.spacing(i)/ls.spacing(i - 1) > .99);
  }
}
