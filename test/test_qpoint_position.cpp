#include <catch2/catch.hpp>
#include <qpoint_position.hpp>
#include <Equidistant.hpp>
#include <Deformed_element.hpp>
#include <cartdgConfig.hpp>

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
  REQUIRE(qpoint_position(elem, basis, row_size*(row_size - 1) + row_size/2)[0] == Approx(1. - 0.2/2).scale(1.));
  REQUIRE(qpoint_position(elem, basis, row_size*(row_size - 1) + row_size/2)[1] == Approx(0.7/2).scale(1.));
  // test the qpoint at the middle of the element (the mean of the vertex positions)
  REQUIRE(qpoint_position(elem, basis, params.n_qpoint()/2)[0] == Approx(.5 - 0.2/4).scale(1.));
  REQUIRE(qpoint_position(elem, basis, params.n_qpoint()/2)[1] == Approx(.5 - 0.3/4).scale(1.));
}

TEST_CASE("jacobian")
{
  const int row_size = 3;
  static_assert (row_size <= cartdg::config::max_row_size);
  cartdg::Equidistant basis {row_size};
  double* jac;
  cartdg::Storage_params params2 {2, 4, 2, row_size};
  cartdg::Deformed_element elem0 {params2, {0, 0}};
  cartdg::Deformed_element elem1 {params2, {1, 1}};
  elem0.vertex(3).pos = {0.8*0.2, 0.8*0.2, 0.};
  elem1.node_adjustments()[6 + 1] = 0.1;
  // jacobian is correct
  set_jacobian(elem0, basis);
  jac = elem0.jacobian();
  REQUIRE(jac[0*9    ] == Approx(1.));
  REQUIRE(jac[1*9    ] == Approx(0.));
  REQUIRE(jac[2*9    ] == Approx(0.));
  REQUIRE(jac[3*9    ] == Approx(1.));
  REQUIRE(jac[0*9 + 6] == Approx(1.));
  REQUIRE(jac[1*9 + 6] == Approx(-0.2));
  REQUIRE(jac[2*9 + 6] == Approx(0.));
  REQUIRE(jac[3*9 + 6] == Approx(0.8));
  REQUIRE(jac[0*9 + 8] == Approx(0.8));
  REQUIRE(jac[1*9 + 8] == Approx(-0.2));
  REQUIRE(jac[2*9 + 8] == Approx(-0.2));
  REQUIRE(jac[3*9 + 8] == Approx(0.8));
  set_jacobian(elem1, basis);
  jac = elem1.jacobian();
  REQUIRE(jac[0*9 + 5] == Approx(1.));
  REQUIRE(jac[1*9 + 5] == Approx(0.));
  REQUIRE(jac[2*9 + 5] == Approx(0.));
  REQUIRE(jac[3*9 + 5] == Approx(0.9));
  // time step is correct
  // At corner 3, diagonal is locally scaled by 1 - 2*(1 - 0.8) = 0.6 = (min singular value)
  REQUIRE(elem0.vertex_time_step_scale()[0] == 1.);
  REQUIRE(elem0.vertex_time_step_scale()[3] == Approx(0.6));

  cartdg::Storage_params params3 {2, 5, 3, row_size};
  cartdg::Deformed_element elem2 {params3, {0, 0, 0}};
  elem2.vertex(7).pos = {0.8*0.2, 0.8*0.2, 0.8*0.2};
  set_jacobian(elem2, basis);
  jac = elem2.jacobian();
  REQUIRE(jac[0] == 1.);
  REQUIRE(jac[0*27 + 26] == Approx( 0.8));
  REQUIRE(jac[1*27 + 26] == Approx(-0.2));
  REQUIRE(jac[2*27 + 26] == Approx(-0.2));
  REQUIRE(jac[7*27 + 26] == Approx(-0.2));
  REQUIRE(jac[8*27 + 26] == Approx( 0.8));
}
