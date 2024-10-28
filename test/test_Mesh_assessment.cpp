#include <catch2/catch_all.hpp>
#include <hexed/Mesh_assessment.hpp>
#include <hexed/Array.hpp>

TEST_CASE("Mesh_assessment") {
  SECTION("2D") {
    hexed::Array<double> vert_arr({2, 2, 3});
    vert_arr = 0;
    vert_arr(1)(0)[0] = .25;
    vert_arr(1)(1)[0] = .25;
    vert_arr(0)(1)[1] = .25;
    vert_arr(1)(1)[1] = .375;
    hexed::Array<double> flat_arr = vert_arr.reshaped({hexed::whatever, 3});
    hexed::next::Sequence<hexed::Mat<3>> vert_seq {
      [&](hexed::Int i)->hexed::Mat<3> {return flat_arr(i).vector();},
      []()->hexed::Int {return 4;},
    };
    hexed::Mesh_assessment ma;
    ma = hexed::Mesh_assessment(vert_seq, 0, 0);
    REQUIRE(ma.orthogonality == Catch::Approx(1.));
    REQUIRE_THAT(ma.edge_lengths, Catch::Matchers::RangeEquals(hexed::Mat<3>{.25, .25, 0.},
                                                               hexed::math::Approx_equal(0, 1e-10)));
    ma = hexed::Mesh_assessment(vert_seq, 1, 0);
    REQUIRE(ma.orthogonality == Catch::Approx(2./std::sqrt(5.)));
    REQUIRE_THAT(ma.edge_lengths, Catch::Matchers::RangeEquals(hexed::Mat<3>{.125*std::sqrt(5.), .25, 0.},
                                                               hexed::math::Approx_equal(0, 1e-10)));
  }
  SECTION("3D") {
    hexed::Array<double> vert_arr({2, 2, 2, 3});
    vert_arr = 0;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          vert_arr(i)(j)(k) = hexed::Array<double>::make(i, j, k);
        }
      }
    }
    vert_arr(1)(0)(0)[0] = .8;
    hexed::Array<double> flat_arr = vert_arr.reshaped({hexed::whatever, 3});
    hexed::next::Sequence<hexed::Mat<3>> vert_seq {
      [&](hexed::Int i)->hexed::Mat<3> {return flat_arr(i).vector();},
      []()->hexed::Int {return 8;},
    };
    hexed::Mesh_assessment ma;
    ma = hexed::Mesh_assessment(vert_seq, 4, 4);
    REQUIRE(ma.orthogonality == Catch::Approx(1/(1 + .2*.2)));
    REQUIRE_THAT(ma.edge_lengths,
                 Catch::Matchers::RangeEquals(hexed::Mat<3>{.8, std::sqrt(1 + .2*.2), std::sqrt(1 + .2*.2)},
                                              hexed::math::Approx_equal(0, 1e-10)));
    ma = hexed::Mesh_assessment(vert_seq, 6, 7);
    REQUIRE(ma.orthogonality == Catch::Approx(1/std::sqrt(1 + .2*.2)));
    REQUIRE_THAT(ma.edge_lengths,
                 Catch::Matchers::RangeEquals(hexed::Mat<3>{1., std::sqrt(1 + .2*.2), 1.},
                                              hexed::math::Approx_equal(0, 1e-10)));
  }
}
