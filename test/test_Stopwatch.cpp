#include <thread>
#include <catch2/catch.hpp>
#include <hexed/Stopwatch.hpp>

TEST_CASE("Stopwatch")
{
  hexed::Stopwatch watch;
  REQUIRE(watch.n_calls() == 0);
  REQUIRE(watch.time() == 0.);
  REQUIRE(!watch.running());
  watch.start();
  REQUIRE_THROWS(watch.start()); // can't start a pausewatch that's already running
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  REQUIRE(watch.running());
  REQUIRE(watch.n_calls() == 0);
  REQUIRE(watch.time() == Approx(0).margin(1e-3));
  watch.pause();
  REQUIRE(!watch.running());
  REQUIRE_THROWS(watch.pause()); // can only pause a running pausewatch
  REQUIRE(watch.n_calls() == 1);
  REQUIRE(watch.n_calls() == 1); // n_calls is incremented by `start` not `n_calls`
  REQUIRE(watch.time() == Approx(1e-2).epsilon(1e-1));
  watch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(20));
  watch.pause();
  REQUIRE(watch.time() == Approx(3e-2).epsilon(1e-1));
  REQUIRE(watch.n_calls() == 2);
  watch.reset();
  REQUIRE(watch.time() == 0.);
  REQUIRE(watch.n_calls() == 0);
}
