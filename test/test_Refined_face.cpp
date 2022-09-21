#include <catch2/catch.hpp>

#include <hexed/Storage_params.hpp>
#include <hexed/Refined_face.hpp>

TEST_CASE("Refined_face")
{
  double coarse {0.}; // usually this should be an array, but one double is fine for testing
  double* coarse_start {&coarse};
  hexed::Storage_params params {3, 5, 3, 2};
  hexed::Refined_face ref_face {params, coarse_start};
  for (int i_face = 0; i_face < 4; ++i_face) {
    ref_face.fine_face(i_face)[0] = i_face;
    ref_face.fine_face(i_face)[5*3] = i_face;
  }
  for (int i_face = 0; i_face < 4; ++i_face) REQUIRE(ref_face.fine_face(i_face)[0] == i_face);
  REQUIRE(ref_face.coarse_face() == coarse_start);
}
