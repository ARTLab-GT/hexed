#include <catch2/catch.hpp>
#include <Element_container.hpp>

TEST_CASE("Specific_container<Deformed_element>")
{
  int row_size = 2;
  cartdg::Storage_params params {3, 5, 3, row_size};
  cartdg::Complete_element_container<cartdg::Deformed_element> ctn {params, 0.3};
  // test that elements are constructed correctly
  int sn0 = ctn.emplace(2, {2, 3});
  cartdg::Deformed_element& elem = ctn.at(2, sn0);
  REQUIRE(elem.vertex(0).pos[1] == Approx(0.3*3./4.));
  REQUIRE(elem.vertex(7).pos[1] == Approx(0.3*4./4.));
}
