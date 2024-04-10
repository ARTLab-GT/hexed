#include <catch2/catch_all.hpp>
#include <hexed/mesh_objects.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Element_new")
{
  constexpr int rs = hexed::config::max_row_size;
  hexed::Gauss_legendre basis(rs);
  hexed::Array<int> pos_ind({2}); pos_ind[0] = 1; pos_ind[1] = 3;
  hexed::Array<double> og({2}); og[0] = .1; og[1] = .2;
  SECTION("invalid construction") {
    REQUIRE_THROWS(hexed::Element_new({2, 4, 2, rs - 1}, false, 2, pos_ind, .8, og, basis));
  }
  SECTION("Cartesian element") {
    hexed::Element_new elem({2, 4, 2, rs}, false, 2, pos_ind, .8, og, basis);
    REQUIRE(elem.deformed() == false);
    REQUIRE(elem.ref_level() == 2);
    REQUIRE(elem.aniso_ref_level() == 0);
    REQUIRE(elem.storage_params().n_dim == 2);
    REQUIRE(elem.position_index().size() == 2);
    REQUIRE(elem.position_index()[1] == 3);
    REQUIRE(elem.root_size() == Catch::Approx(.8));
    REQUIRE(elem.nominal_size() == Catch::Approx(.2));
    REQUIRE(elem.nominal_position()[0] == Catch::Approx(.3));
    REQUIRE(elem.nominal_position()[1] == Catch::Approx(.8));
    REQUIRE(elem.basis().row_size == rs);
  }
}
