#include <catch2/catch_all.hpp>
#include <hexed/Mesh_assessment.hpp>
#include <hexed/Array.hpp>

const double grad_scale = 1e-6;

#define TEST_GRADIENT \
  for (int i = 0; i < vert_seq.size(); ++i) { \
    for (int j = 0; j < vert_seq.size(); ++j) { \
      ma = hexed::Mesh_assessment(vert_seq, i, j); \
      hexed::Mat<3> orth = ma.orthogonality; \
      hexed::Mat<3> lengths = ma.edge_lengths; \
      hexed::Array<double> old_pos {flat_arr(j).copy()}; \
      hexed::Mat<3> vec; \
      for (int i = 0; i < 3; ++i) vec(i) = (rand()%2000 - 1000)*1e-3; \
      flat_arr(j).vector() += grad_scale*vec; \
      ma = hexed::Mesh_assessment(vert_seq, i, j); \
      flat_arr(j) = old_pos; \
      for (int i_dim = 0; i_dim < 3; ++i_dim) { \
        REQUIRE((ma.orthogonality(i_dim) - orth(i_dim))/grad_scale \
                == Catch::Approx(vec.dot(ma.grad_orth(i_dim, hexed::all))).margin(1e-4)); \
        REQUIRE((ma.edge_lengths(i_dim) - lengths(i_dim))/grad_scale \
                == Catch::Approx(vec.dot(ma.grad_lengths(i_dim, hexed::all))).margin(1e-4)); \
      } \
    } \
  } \

TEST_CASE("Mesh_assessment") {
  srand(1011);
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
    REQUIRE_THAT(ma.orthogonality, Catch::Matchers::RangeEquals(hexed::Mat<3>{1., 1., 1.},
                                                                hexed::math::Approx_equal(0, 1e-10)));
    REQUIRE_THAT(ma.edge_lengths, Catch::Matchers::RangeEquals(hexed::Mat<3>{.25, .25, 0.},
                                                               hexed::math::Approx_equal(0, 1e-10)));
    ma = hexed::Mesh_assessment(vert_seq, 1, 0);
    REQUIRE_THAT(ma.orthogonality, Catch::Matchers::RangeEquals(hexed::Mat<3>{2./std::sqrt(5.), 2./std::sqrt(5.), 1.},
                                                                hexed::math::Approx_equal(0, 1e-10)));
    REQUIRE_THAT(ma.edge_lengths, Catch::Matchers::RangeEquals(hexed::Mat<3>{.125*std::sqrt(5.), .25, 0.},
                                                               hexed::math::Approx_equal(0, 1e-10)));
    TEST_GRADIENT
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
    REQUIRE_THAT(ma.edge_lengths,
                 Catch::Matchers::RangeEquals(hexed::Mat<3>{.8, std::sqrt(1 + .2*.2), std::sqrt(1 + .2*.2)},
                                              hexed::math::Approx_equal(0, 1e-10)));
    ma = hexed::Mesh_assessment(vert_seq, 6, 7);
    REQUIRE_THAT(ma.edge_lengths,
                 Catch::Matchers::RangeEquals(hexed::Mat<3>{1., std::sqrt(1 + .2*.2), 1.},
                                              hexed::math::Approx_equal(0, 1e-10)));
    TEST_GRADIENT
    vert_arr(1)(1)(0)[0] = .8;
    ma = hexed::Mesh_assessment(vert_seq, 4, 4);
    double orth = 1./std::sqrt(1. + .2*.2);
    REQUIRE_THAT(ma.orthogonality,
                 Catch::Matchers::RangeEquals(hexed::Mat<3>{orth, 1., orth},
                                              hexed::math::Approx_equal(0, 1e-10)));
  }
}
