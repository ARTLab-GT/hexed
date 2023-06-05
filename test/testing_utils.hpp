#include <catch2/catch_all.hpp>
#include <hexed/Vector_view.hpp>
#include <hexed/Deformed_element.hpp>

typedef hexed::Vector_view<hexed::Element&, std::unique_ptr<hexed::Element>, &hexed::ptr_convert<hexed::Element&, std::unique_ptr<hexed::Element>>> car_elem_view;
typedef hexed::Vector_view<hexed::Deformed_element&, std::unique_ptr<hexed::Deformed_element>, &hexed::ptr_convert<hexed::Deformed_element&, std::unique_ptr<hexed::Deformed_element>>> def_elem_view;

// useful for a few very specific tests involving `hexed::Vertex`s
inline void assert_equal(std::array<double, 3> computed, std::array<double, 3> correct)
{
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    CHECK(computed[i_dim] == Catch::Approx(correct[i_dim]).margin(1e-14));
  }
}

template <typename T, typename U>
void require_sequence_equal(T sequence0, U sequence1)
{
  REQUIRE(sequence0.size() == sequence1.size());
  for (int i = 0; i < sequence0.size(); ++i) {
    CHECK(sequence0[i] == Catch::Approx(sequence1[i]));
  }
}
