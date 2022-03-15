#include <catch2/catch.hpp>
#include <connection.hpp>

TEST_CASE("Cartesian_element_connection")
{
  cartdg::Storage_params params {3, 5, 3, 6};
  cartdg::Element elem0 (params);
  cartdg::Element elem1 (params);
  cartdg::Cartesian_element_connection con ({&elem0, &elem1}, 1);
  REQUIRE(con.i_dim() == 1);
  REQUIRE(&con.element(0) == &elem0);
  REQUIRE(&con.element(1) == &elem1);
  REQUIRE(con.face(0) == elem0.face() + 3*5*36);
  REQUIRE(con.face(1) == elem1.face() + 1*5*36);
}

TEST_CASE("Deformed_element_connection")
{
  cartdg::Storage_params params {3, 5, 3, 6};
  cartdg::Deformed_element elem0 (params);
  cartdg::Deformed_element elem1 (params);
  cartdg::Deformed_element_connection con ({&elem0, &elem1}, {2, 1}, {0, 1});
  REQUIRE(con.i_dim(0) == 2);
  REQUIRE(con.i_dim(1) == 1);
  REQUIRE(con.face_sign(0) == 0);
  REQUIRE(con.face_sign(1) == 1);
  REQUIRE(&con.element(0) == &elem0);
  REQUIRE(&con.element(1) == &elem1);
  con.jacobian()[3*3*6*6 - 1] = 1; // check that jacobian is big enough (otherwise segfault)
  REQUIRE(con.face(0) == elem0.face() + 4*5*36);
  REQUIRE(con.face(1) == elem1.face() + 3*5*36);
}
