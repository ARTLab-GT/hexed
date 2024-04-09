#include <catch2/catch_all.hpp>
#include <hexed/mesh_objects.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Element_new")
{
  constexpr int rs = hexed::config::max_row_size;
  hexed::Gauss_legendre basis(rs);
  hexed::Array<int> pos_ind({2}); pos_ind[0] = 1; pos_ind[1] = 3;
  hexed::Array<double> og({2}); og[0] = .1; og[1] = .2;
  hexed::Element_new elem({2, 4, 2, rs}, false, 2, pos_ind, .8, og, basis);
  REQUIRE_THROWS(hexed::Element_new({2, 4, 2, rs - 1}, false, 2, pos_ind, .8, og, basis));
}
