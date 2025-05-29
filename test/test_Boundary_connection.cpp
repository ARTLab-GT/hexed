#include <catch2/catch_all.hpp>
#include <hexed/Boundary_connection.hpp>

TEST_CASE("Boundary_connection") {
  hexed::Storage_params params {2, 5, 3, 4};
  hexed::Face f0(params, 2, 0, false);
  hexed::Face f1(params, 1, 1, true);
  hexed::next::Boundary_connection bc0(f0, 5, 1);
  REQUIRE(bc0.boundary_condition() == 5);
  REQUIRE(&bc0.neighbor_connection().face(0) == &f0);
  REQUIRE(&bc0.inside() == &f0);
  REQUIRE(bc0.neighbor_connection().face(1).storage_params() == f0.storage_params());
  REQUIRE(bc0.neighbor_connection().face(1).i_dim() == 2);
  REQUIRE(bc0.neighbor_connection().face(1).sign() == 1);
  REQUIRE(bc0.neighbor_connection().face(1).is_deformed() == false);
  REQUIRE(bc0.ghost().associated());
  REQUIRE(bc0.ghost().boundary_connection() == &bc0);
  REQUIRE(f0.connected());
  REQUIRE_THROWS(hexed::next::Boundary_connection(f0, 5, 1));
  hexed::next::Boundary_connection bc1(f1, 0, 6);
  REQUIRE(bc1.ghost().i_dim() == 1);
  REQUIRE(bc1.ghost().sign() == 0);
  REQUIRE(bc1.ghost().is_deformed() == true);
  REQUIRE_THAT(bc0.normal().shape(), Catch::Matchers::RangeEquals(std::vector<hexed::Int>{3, 16}));
  REQUIRE_THAT(bc1.normal().shape(), Catch::Matchers::RangeEquals(std::vector<hexed::Int>{3, 16}));
  f1.normal()(0) = .1;
  f1.normal()(1) = .2;
  f1.normal()(2) = .3;
  for (int i = 0; i < 16; ++i) {
    REQUIRE(bc0.normal()(0)[i] == Catch::Approx(0.).scale(1.));
    REQUIRE(bc0.normal()(1)[i] == Catch::Approx(0.).scale(1.));
    REQUIRE(bc0.normal()(2)[i] == Catch::Approx(1.).scale(1.));
    REQUIRE(bc1.normal()(0)[i] == Catch::Approx(.1));
    REQUIRE(bc1.normal()(1)[i] == Catch::Approx(.2));
    REQUIRE(bc1.normal()(2)[i] == Catch::Approx(.3));
  }
}
