#include <thread>
#include <catch2/catch_all.hpp>
#include <hexed/Stopwatch.hpp>

TEST_CASE("Stopwatch")
{
  hexed::Stopwatch watch;
  REQUIRE(watch.n_calls() == 0);
  REQUIRE(watch.time() == 0.);
  REQUIRE(!watch.running());
  watch.start();
  REQUIRE_THROWS(watch.start()); // can't start a pausewatch that's already running
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  REQUIRE(watch.running());
  REQUIRE(watch.n_calls() == 0);
  REQUIRE(watch.time() == Catch::Approx(0).margin(.01));
  watch.pause();
  REQUIRE(!watch.running());
  REQUIRE_THROWS(watch.pause()); // can only pause a running pausewatch
  REQUIRE(watch.n_calls() == 1);
  REQUIRE(watch.n_calls() == 1); // n_calls is incremented by `start` not `n_calls`
  REQUIRE(watch.time() == Catch::Approx(.1).epsilon(1e-1));
  watch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  watch.pause();
  REQUIRE(watch.time() == Catch::Approx(.3).epsilon(1e-1));
  REQUIRE(watch.n_calls() == 2);
  watch.reset();
  REQUIRE(watch.time() == 0.);
  REQUIRE(watch.n_calls() == 0);
  hexed::Stopwatch watch1;
  watch1.start();
  watch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  REQUIRE_THROWS(watch + watch1);
  watch1.pause();
  watch.pause();
  hexed::Stopwatch sum = watch + watch1;
  REQUIRE(sum.time() == Catch::Approx(.4).epsilon(.1));
  REQUIRE(sum.n_calls() == 2);
}
