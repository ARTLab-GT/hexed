#include <catch2/catch.hpp>
#include <Element_container.hpp>

TEST_CASE("Specific_container<Deformed_element>")
{
  int row_size = 2;
  cartdg::Storage_params params {3, 5, 3, row_size};
  cartdg::Complete_element_container<cartdg::Deformed_element> ctn {params, 0.3};
  // test access by refinement level and serial number
  int sn0 = ctn.emplace(2, {2, 3});
  int sn1 = ctn.emplace(2, {2, 1});
  int sn2 = ctn.emplace(3, {2, 1});
  int sn3 = ctn.emplace(2, {2, 1});
  REQUIRE(&ctn.at(2, sn0) != &ctn.at(2, sn1));
  REQUIRE(&ctn.at(2, sn3) != &ctn.at(2, sn1));
  REQUIRE_THROWS(ctn.at(2, std::abs(sn0) + std::abs(sn1) + 1));
  REQUIRE_THROWS(ctn.at(3, sn0));
  REQUIRE(&ctn.at(3, sn2) != &ctn.at(2, sn1));
  // test that elements are constructed correctly
  cartdg::Deformed_element& elem = ctn.at(2, sn0);
  REQUIRE(elem.vertex(0).pos[1] == Approx(0.3*3./4.));
  REQUIRE(elem.vertex(7).pos[1] == Approx(0.3*4./4.));
  REQUIRE(ctn.at(3, sn2).vertex(0).pos[0] == Approx(0.3*2./8.));
  // test that `ctn` purports to contain 4 `Deformed_element`s
  REQUIRE(ctn.elements().size() == 4);
  // test that each of the `Deformed_element` constructed above appears exactly once in the vector view
  cartdg::Deformed_element* ptrs [] {&ctn.at(2, sn0), &ctn.at(2, sn1), &ctn.at(3, sn2), &ctn.at(2, sn3)};
  for (int i_elem = 0; i_elem < 4; ++i_elem) {
    int count = 0;
    for (int j_elem = 0; j_elem < 4; ++j_elem) {
      cartdg::Deformed_element& elem {ctn.elements()[j_elem]};
      if (&elem == ptrs[i_elem]) ++count;
    }
    REQUIRE(count == 1);
  }
}