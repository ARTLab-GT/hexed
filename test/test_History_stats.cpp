#include <random>
#include <catch2/catch_all.hpp>
#include <hexed/History_stats.hpp>

TEST_CASE("History_stats") {
  hexed::History_stats stats(.001);
  REQUIRE(stats.n_sample() == 0);
  REQUIRE(stats.last_iter() == -1);
  REQUIRE(std::isnan(stats.last_value()));
  REQUIRE(std::isnan(stats.mean()));
  REQUIRE(std::isnan(stats.std_dev()));
  REQUIRE(std::isnan(stats.deriv()));
  REQUIRE(std::isnan(stats.deriv_std_dev()));

  SECTION("first iterations") {
    stats.add_sample(2, .1);
    REQUIRE(stats.n_sample() == 1);
    REQUIRE(stats.last_value() == Catch::Approx(.1));
    REQUIRE(stats.mean() == Catch::Approx(.1));
    REQUIRE(std::isnan(stats.std_dev()));
    REQUIRE(std::isnan(stats.deriv()));
    REQUIRE(std::isnan(stats.deriv_std_dev()));
    stats.add_sample(10, .2);
    REQUIRE(stats.n_sample() == 2);
    REQUIRE(stats.std_dev() == Catch::Approx(0.).scale(1.));
    REQUIRE(stats.deriv() == Catch::Approx(.1/8));
    REQUIRE(std::isnan(stats.deriv_std_dev()));
    stats.add_sample(20, .3);
    REQUIRE(stats.n_sample() == 3);
    REQUIRE(std::isfinite(stats.deriv_std_dev()));
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution normal(0.);
  hexed::Int n_iter = 10'000'000;
  for (hexed::Int i = 0; i <= n_iter; i += 10) {
    stats.add_sample(i, 1 - std::exp(-double(i)/n_iter) + 1e-2/(1 + double(i)/n_iter)*normal(gen));
  }
  REQUIRE(stats.n_sample() == n_iter/10 + 1);
  double curr_mean = 1 - std::exp(-1.);
  REQUIRE(stats.last_iter() == n_iter);
  REQUIRE(stats.last_value() == Catch::Approx(curr_mean).epsilon(.05));
  REQUIRE(stats.mean() == Catch::Approx(curr_mean).epsilon(.05));
  REQUIRE(stats.std_dev() == Catch::Approx(1e-2/2.).epsilon(.05));
  REQUIRE(stats.deriv()*n_iter == Catch::Approx(std::exp(-1.)).epsilon(.05));
  REQUIRE(stats.deriv_std_dev() == Catch::Approx(-1e-2/4*1./n_iter).epsilon(.05));
}
