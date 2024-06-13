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

    REQUIRE_THAT(elem.full_state().shape()        , Catch::Matchers::RangeEquals(std::vector<int>{11 + 2*rs, rs, rs}));
    REQUIRE_THAT(elem.flow_state().shape()        , Catch::Matchers::RangeEquals(std::vector<int>{4, rs, rs}));
    REQUIRE_THAT(elem.ltss().shape()              , Catch::Matchers::RangeEquals(std::vector<int>{rs, rs}));
    REQUIRE_THAT(elem.bulk_art_visc().shape()     , Catch::Matchers::RangeEquals(std::vector<int>{rs, rs}));
    REQUIRE_THAT(elem.laplacian_art_visc().shape(), Catch::Matchers::RangeEquals(std::vector<int>{rs, rs}));
    REQUIRE_THAT(elem.art_visc_forcing().shape()  , Catch::Matchers::RangeEquals(std::vector<int>{4, rs, rs}));
    REQUIRE_THAT(elem.advection_state().shape()   , Catch::Matchers::RangeEquals(std::vector<int>{rs, rs, rs}));
    REQUIRE_THAT(elem.cache().shape()             , Catch::Matchers::RangeEquals(std::vector<int>{std::max(4, rs), rs, rs}));
    REQUIRE(elem.flow_state().data()         == elem.full_state().data());
    REQUIRE(elem.ltss().data()               == elem.full_state()(4).data());
    REQUIRE(elem.bulk_art_visc().data()      == elem.full_state()(5).data());
    REQUIRE(elem.laplacian_art_visc().data() == elem.full_state()(6).data());
    REQUIRE(elem.art_visc_forcing().data()   == elem.full_state()(7).data());
    REQUIRE(elem.advection_state().data()    == elem.full_state()(11).data());
    REQUIRE(elem.cache().data()              == elem.full_state()(11 + rs).data());
  }
}
