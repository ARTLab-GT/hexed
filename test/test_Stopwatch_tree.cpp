#include <thread>
#include <catch2/catch.hpp>
#include <Stopwatch_tree.hpp>

TEST_CASE("Stopwatch_tree")
{
  cartdg::Stopwatch_tree tree ("unit0");
  {
    cartdg::Stopwatch_tree::Measurement mmt (tree, 3);
    REQUIRE_THROWS(cartdg::Stopwatch_tree::Measurement (tree, 3));
    REQUIRE(tree.work_units_completed() == 0);
    REQUIRE(tree.time() == Approx(0).margin(1e-3));
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }
  REQUIRE(tree.work_units_completed() == 3);
  REQUIRE(tree.time() == Approx(2e-2).epsilon(1e-1));
}
