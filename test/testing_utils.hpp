#include <catch2/catch.hpp>

inline void assert_equal(std::array<double, 3> computed, std::array<double, 3> correct)
{
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    CHECK(computed[i_dim] == Approx(correct[i_dim]).margin(1e-14));
  }
}
