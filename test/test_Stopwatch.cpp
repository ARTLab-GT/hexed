#include <catch2/catch.hpp>
#include <Stopwatch.hpp>

TEST_CASE("Stopwatch")
{
  cartdg::Stopwatch watch;
  watch.start();
  REQUIRE_THROWS(watch.start()); // can't start a pausewatch that's already running
  watch.pause();
  REQUIRE_THROWS(watch.pause()); // can only pause a running pausewatch
  watch.start();
  watch.pause();
}
