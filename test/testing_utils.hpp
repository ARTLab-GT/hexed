#include <catch2/catch.hpp>
#include <Vector_view.hpp>
#include <Deformed_element.hpp>

typedef cartdg::Vector_view<cartdg::Element&, std::unique_ptr<cartdg::Element>, &cartdg::ptr_convert<cartdg::Element&, std::unique_ptr<cartdg::Element>>> car_elem_view;
typedef cartdg::Vector_view<cartdg::Deformed_element&, std::unique_ptr<cartdg::Deformed_element>, &cartdg::ptr_convert<cartdg::Deformed_element&, std::unique_ptr<cartdg::Deformed_element>>> def_elem_view;

// useful for a few very specific tests involving `cartdg::Vertex`s
inline void assert_equal(std::array<double, 3> computed, std::array<double, 3> correct)
{
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    CHECK(computed[i_dim] == Approx(correct[i_dim]).margin(1e-14));
  }
}
