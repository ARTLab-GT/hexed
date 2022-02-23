#include <catch2/catch.hpp>
#include <Element_container.hpp>

TEST_CASE("Specific_container<Deformed_element>")
{
  int row_size = 2;
  cartdg::Storage_params params {3, 5, 3, row_size};
  cartdg::Specific_container<cartdg::Deformed_element> cont {params, 0.3};
}
