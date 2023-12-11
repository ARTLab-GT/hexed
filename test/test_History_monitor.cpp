#include <catch2/catch_all.hpp>
#include <hexed/History_monitor.hpp>

TEST_CASE("History_monitor")
{
  hexed::History_monitor hist(.3, 10);
  REQUIRE(hist.max() ==  std::numeric_limits<double>::max());
  REQUIRE(hist.min() == -std::numeric_limits<double>::max());
  hist.add_sample(0, 1.);
  REQUIRE(hist.max() == Catch::Approx(1.));
  REQUIRE(hist.min() == Catch::Approx(1.));
  for (int i = 1; i < 10; ++i) hist.add_sample(i, 1./(i + 1));
  CHECK(hist.min() == Catch::Approx(.1));
  CHECK(hist.max() == Catch::Approx(1./8.));
  for (int i = 10; i < 1000; ++i) hist.add_sample(i, 1./(i + 1));
  REQUIRE(hist.min() == Catch::Approx(1e-3).margin(1e-6*(.3/10)*1000));
  REQUIRE(hist.max() == Catch::Approx(1/.7e3).margin(1/.7e3/.7e3*(.3/10)*1000));
}
