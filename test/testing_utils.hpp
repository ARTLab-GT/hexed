#include <catch2/catch.hpp>
#include <Vector_view.hpp>
#include <Deformed_element.hpp>

typedef hexed::Vector_view<hexed::Element&, std::unique_ptr<hexed::Element>, &hexed::ptr_convert<hexed::Element&, std::unique_ptr<hexed::Element>>> car_elem_view;
typedef hexed::Vector_view<hexed::Deformed_element&, std::unique_ptr<hexed::Deformed_element>, &hexed::ptr_convert<hexed::Deformed_element&, std::unique_ptr<hexed::Deformed_element>>> def_elem_view;

// useful for a few very specific tests involving `hexed::Vertex`s
inline void assert_equal(std::array<double, 3> computed, std::array<double, 3> correct)
{
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    CHECK(computed[i_dim] == Approx(correct[i_dim]).margin(1e-14));
  }
}
