#include <catch2/catch_all.hpp>
#include <hexed/History_monitor.hpp>

TEST_CASE("History_monitor") {
  hexed::History_monitor hist(.3, 10);
  for (int i = 0; i < 2; ++i) {
    REQUIRE(hist.max() ==  std::numeric_limits<double>::max());
    REQUIRE(hist.min() == -std::numeric_limits<double>::max());
    hist.add_sample(0, 1.);
    REQUIRE(hist.max() == Catch::Approx(1.));
    REQUIRE(hist.min() == Catch::Approx(1.));
    for (int j = 1; j < 10; ++j) hist.add_sample(j, 1./(j + 1));
    REQUIRE(hist.min() == Catch::Approx(.1));
    REQUIRE(hist.max() == Catch::Approx(1./8.));
    for (int j = 10; j < 1000; ++j) hist.add_sample(j, 1./(j + 1));
    REQUIRE(hist.min() == Catch::Approx(1e-3).margin(1e-6*(.3/10)*1000));
    REQUIRE(hist.max() == Catch::Approx(1/.7e3).margin(1/.7e3/.7e3*(.3/10)*1000));
    hist.clear();
  }
}
