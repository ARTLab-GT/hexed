#include <catch2/catch.hpp>
#include <Solution.hpp>

// useful for a few very specific tests involving `cartdg::Vertex`s
inline void assert_equal(std::array<double, 3> computed, std::array<double, 3> correct)
{
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    CHECK(computed[i_dim] == Approx(correct[i_dim]).margin(1e-14));
  }
}

/*
 * Sets up a grid for performing integrated tests with deformed
 * and regular elements. It is designed to incur as many of the
 * programming challenges associated with deformed grids as possible
 * to verify proper treatment. See implementation for additional details.
 * Returns a `std::uniqe_ptr` to avoid dealing with the issue of copy/move
 * semantics for `cartdg::Solution`.
 */
std::unique_ptr<cartdg::Solution> deformed_test_setup_2d(bool perturb_face = true);
