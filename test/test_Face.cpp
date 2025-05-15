#include <catch2/catch_all.hpp>
#include <hexed/Face.hpp>

TEST_CASE("Face") {
  hexed::Storage_params params {2, 5, 3, 2};
  hexed::Face f(params, 0, 1);
  REQUIRE(f.i_dim() == 0);
  REQUIRE(f.sign() == 1);
  REQUIRE(f.storage_params().row_size == 2);
  REQUIRE(!f.connected());
  REQUIRE(!f.associated());
  REQUIRE(f.element() == nullptr);
  hexed::Element elem0(params);
  hexed::Element elem1(params);
  f.associate(elem0);
  REQUIRE(f.element() == &elem0);
  REQUIRE(f.associated());
  REQUIRE_THROWS(f.associate(elem1));
}
