#include <catch2/catch.hpp>
#include <hexed/Layer_sequence.hpp>

TEST_CASE("Layer_sequence")
{
  hexed::Layer_sequence ls(1e-6, 4);
  REQUIRE(ls.cumulative_height(ls.n_layers()) == Approx(1.));
  for (int i = 0; i < 4; ++i) REQUIRE(ls.spacing(i) == Approx(1e-6));
  REQUIRE(ls.cumulative_height(0) == Approx(0.).scale(1e-6));
  REQUIRE(ls.cumulative_height(4) == Approx(0.).scale(4e-6));
  REQUIRE(ls.spacing(ls.n_layers() - 1) > .5);
  for (int i = 1; i < ls.n_layers(); ++i) {
    REQUIRE(ls.spacing(i)/ls.spacing(i - 1) < 2);
  }
}
