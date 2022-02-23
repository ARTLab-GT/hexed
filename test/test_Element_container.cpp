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
  REQUIRE(&ctn.at(2, sn0).element != &ctn.at(2, sn1).element);
  REQUIRE(&ctn.at(2, sn3).element != &ctn.at(2, sn1).element);
  REQUIRE_THROWS(ctn.at(2, std::abs(sn0) + std::abs(sn1) + 1));
  REQUIRE_THROWS(ctn.at(3, sn0));
  REQUIRE(&ctn.at(3, sn2).element != &ctn.at(2, sn1).element);
  // test that elements are constructed correctly
  cartdg::Deformed_element& elem = ctn.at(2, sn0).element;
  REQUIRE(elem.vertex(0).pos[1] == Approx(0.3*3./4.));
  REQUIRE(elem.vertex(7).pos[1] == Approx(0.3*4./4.));
  for (int i = 0; i < 6; ++i) REQUIRE(ctn.at(2, sn0).connectedness[i] == 0);
  REQUIRE(ctn.at(3, sn2).element.vertex(0).pos[0] == Approx(0.3*2./8.));
}
