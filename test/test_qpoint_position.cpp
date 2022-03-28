#include <catch2/catch.hpp>
#include <qpoint_position.hpp>
#include <Equidistant.hpp>
#include <Deformed_element.hpp>

TEST_CASE("qpoint_position")
{
  const int row_size = 5;
  cartdg::Storage_params params {2, 4, 2, row_size};
  cartdg::Deformed_element elem {params};
  elem.vertex(3).pos[0] = 0.8;
  elem.vertex(3).pos[1] = 0.7;
  cartdg::Equidistant basis {row_size};
  // test first and last qpoints
  REQUIRE(qpoint_position(elem, basis, 0)[0] == Approx(0.).scale(1.));
  REQUIRE(qpoint_position(elem, basis, 0)[1] == Approx(0.).scale(1.));
  REQUIRE(qpoint_position(elem, basis, 0)[2] == Approx(0.).scale(1.));
  REQUIRE(qpoint_position(elem, basis, params.n_qpoint() - 1)[0] == Approx(.8).scale(1.));
  REQUIRE(qpoint_position(elem, basis, params.n_qpoint() - 1)[1] == Approx(.7).scale(1.));
  static_assert (row_size%2 == 1); // `row_size` must be odd for the following tests to work
  // test the qpoint at the midpoint of the positive-dimension0 face (the right-hand face)
  REQUIRE(qpoint_position(elem, basis, row_size*(row_size/2) - 1)[0] == Approx(1. - 0.2/2).scale(1.));
  REQUIRE(qpoint_position(elem, basis, row_size*(row_size/2) - 1)[1] == Approx(0.7/2).scale(1.));
  // test the qpoint at the middle of the element (the mean of the vertex positions)
  REQUIRE(qpoint_position(elem, basis, params.n_qpoint()/2)[0] == Approx(.5 - 0.2/4).scale(1.));
  REQUIRE(qpoint_position(elem, basis, params.n_qpoint()/2)[1] == Approx(.5 - 0.3/4).scale(1.));
}
