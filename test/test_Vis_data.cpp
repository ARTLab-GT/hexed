#include <catch2/catch_all.hpp>
#include <hexed/Vis_data.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/config.hpp>

TEST_CASE("Vis_data") {
  const int row_size = hexed::config::max_row_size;
  hexed::Array<double> data({4, row_size, row_size, row_size});
  hexed::Gauss_legendre basis(row_size);
  for (int i = 0; i < row_size; ++i) {
    for (int j = 0; j < row_size; ++j) {
      for (int k = 0; k < row_size; ++k) {
        data(0)(i)(j)[k] = basis.node(i);
        data(1)(i)(j)[k] = basis.node(j);
        data(2)(i)(j)[k] = basis.node(k);
      }
    }
  }
  data(3) = data(0)*data(0) + data(1)*data(1) + data(2)*data(2);
  hexed::Vis_data vis(data, basis);

  auto sample = vis.sample(hexed::Array<double>::make(.1, .4, .7, .2, .2, .6).reshaped({3, 2}));
  REQUIRE_THAT(sample.shape(), Catch::Matchers::RangeEquals(std::vector<int>{4, 2}));
  REQUIRE_THAT(sample.column(0).vector(), Catch::Matchers::RangeEquals(std::vector<double>{.1, .7, .2, .1*.1 + .7*.7 + .2*.2},
                                                                       hexed::math::Approx_equal(1e-6)));
  REQUIRE_THAT(sample.column(1).vector(), Catch::Matchers::RangeEquals(std::vector<double>{.4, .2, .6, .4*.4 + .2*.2 + .6*.6},
                                                                       hexed::math::Approx_equal(1e-6)));

  auto interior = vis.interior(21);
  REQUIRE_THAT(interior.shape(), Catch::Matchers::RangeEquals(std::vector<int>{4, 21, 21, 21}));
  REQUIRE(interior(1)(3)(7)[9] == Catch::Approx(7/20.));
  REQUIRE(interior(3)(3)(7)[9] == Catch::Approx((3*3 + 7*7 + 9*9)/20./20.));

  auto edges = vis.edges(21);
  REQUIRE_THAT(edges.shape(), Catch::Matchers::RangeEquals(std::vector<int>{3, 4, 4, 21}));
  REQUIRE(edges(0)(1)(0)[2] == Catch::Approx(.1).scale(1.));
  REQUIRE(edges(0)(1)(1)[2] == Catch::Approx(0.).scale(1.));
  REQUIRE(edges(0)(1)(2)[2] == Catch::Approx(1.).scale(1.));
  REQUIRE(edges(2)(1)(0)[2] == Catch::Approx(0.).scale(1.));
  REQUIRE(edges(2)(1)(1)[2] == Catch::Approx(1.).scale(1.));
  REQUIRE(edges(2)(1)(2)[2] == Catch::Approx(.1).scale(1.));

  auto contour = vis.compute_contour(3, .7, 10, 4, 1e-6);
  auto contour_sample = vis.sample(contour.vert_ref_coords);
  REQUIRE(contour_sample.size());
  for (int i_point = 0; i_point < data.shape()[1]; ++i_point) {
    CHECK(contour_sample(3)[i_point] == Catch::Approx(.7).epsilon(1e-3));
  }
}
