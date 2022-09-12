#include <catch2/catch.hpp>
#include <connection.hpp>

TEST_CASE("direction conversion")
{
  cartdg::Con_dir<cartdg::Element> car {2};
  cartdg::Con_dir<cartdg::Deformed_element> def = car;
  REQUIRE(def.i_dim[0] == 2);
  REQUIRE(def.i_dim[1] == 2);
  REQUIRE(def.face_sign[0] == 1);
  REQUIRE(def.face_sign[1] == 0);
}

TEST_CASE("vertex_inds")
{
  SECTION("Same direction")
  {
    auto inds = cartdg::vertex_inds(3, {{0, 0}, {1, 0}});
    REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 0);
    REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 1);
    REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 2);
    REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 3);
    inds = cartdg::vertex_inds(3, {{1, 1}, {0, 1}});
    REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 2);
    REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 3);
    REQUIRE(inds[0][2] == 4); REQUIRE(inds[1][2] == 6);
    REQUIRE(inds[0][3] == 5); REQUIRE(inds[1][3] == 7);
    inds = cartdg::vertex_inds(3, {{2, 2}, {0, 1}});
    REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 1);
    REQUIRE(inds[0][1] == 2); REQUIRE(inds[1][1] == 3);
    REQUIRE(inds[0][2] == 4); REQUIRE(inds[1][2] == 5);
    REQUIRE(inds[0][3] == 6); REQUIRE(inds[1][3] == 7);
  }

  SECTION("Different direction")
  {
    SECTION("0+ 1+")
    {
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {1, 1}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("0- 1-")
    {
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {0, 0}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 1);
      REQUIRE(inds[0][2] == 2); REQUIRE(inds[1][2] == 4);
      REQUIRE(inds[0][3] == 3); REQUIRE(inds[1][3] == 5);
    }
    SECTION("2+ 1+")
    {
      auto inds = cartdg::vertex_inds(3, {{2, 1}, {1, 1}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("0+ 2+")
    {
      auto inds = cartdg::vertex_inds(3, {{0, 2}, {1, 1}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 1);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 5);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 3);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("2+ 0+")
    {
      auto inds = cartdg::vertex_inds(3, {{2, 0}, {1, 1}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 4);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 6);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 5);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("0+ 1-")
    {
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {1, 0}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 4);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 5);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 1);
    }
    SECTION("0- 1+")
    {
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {0, 1}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 6);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 7);
      REQUIRE(inds[0][2] == 2); REQUIRE(inds[1][2] == 2);
      REQUIRE(inds[0][3] == 3); REQUIRE(inds[1][3] == 3);
    }
    SECTION("1+ 0-")
    {
      auto inds = cartdg::vertex_inds(3, {{1, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 1);
    }
    SECTION("1+ 2-")
    {
      auto inds = cartdg::vertex_inds(3, {{1, 2}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 0);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 4);
    }
    SECTION("2+ 0-")
    {
      auto inds = cartdg::vertex_inds(3, {{2, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 1);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 2);
    }

    SECTION("2D")
    {
      auto inds = cartdg::vertex_inds(2, {{0, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 1);
      inds = cartdg::vertex_inds(2, {{0, 1}, {0, 1}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 3);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 1);
      inds = cartdg::vertex_inds(2, {{1, 0}, {0, 0}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 2); REQUIRE(inds[1][1] == 1);
    }
  }
}

TEST_CASE("Element_face_connection<Element>")
{
  cartdg::Storage_params params {3, 5, 3, 6};
  cartdg::Element elem0 (params);
  cartdg::Element elem1 (params);
  cartdg::Element_face_connection<cartdg::Element> con ({&elem0, &elem1}, cartdg::Con_dir<cartdg::Element>{1});
  REQUIRE(con.direction().i_dim == 1);
  REQUIRE(&con.element(0) == &elem0);
  REQUIRE(&con.element(1) == &elem1);
  REQUIRE(con.face(0) == elem0.face() + 3*5*36);
  REQUIRE(con.face(1) == elem1.face() + 2*5*36);
  REQUIRE(&elem0.vertex(2) == &elem1.vertex(0));
  REQUIRE(&elem0.vertex(7) == &elem1.vertex(5));
  // make sure it didn't just combine all the vertices or something stupid like that
  REQUIRE(&elem0.vertex(0) != &elem1.vertex(2));
}

TEST_CASE("Element_face_connection<Deformed_element>")
{
  cartdg::Storage_params params {3, 5, 3, 6};
  cartdg::Deformed_element elem0 (params);
  cartdg::Deformed_element elem1 (params);
  cartdg::Element_face_connection<cartdg::Deformed_element> con ({&elem0, &elem1}, cartdg::Con_dir<cartdg::Deformed_element>{{2, 1}, {0, 1}});
  REQUIRE(con.direction().i_dim[0] == 2);
  REQUIRE(con.direction().i_dim[1] == 1);
  REQUIRE(con.direction().face_sign[0] == 0);
  REQUIRE(con.direction().face_sign[1] == 1);
  REQUIRE(&con.element(0) == &elem0);
  REQUIRE(&con.element(1) == &elem1);
  con.jacobian()[3*3*6*6 - 1] = 1; // check that jacobian is big enough (otherwise segfault)
  REQUIRE(con.face(0) == elem0.face() + 4*5*36);
  REQUIRE(con.face(1) == elem1.face() + 3*5*36);
  REQUIRE(&elem0.vertex(0) == &elem1.vertex(3));
  REQUIRE(&elem0.vertex(2) == &elem1.vertex(2));
  REQUIRE(&elem0.vertex(4) == &elem1.vertex(7));
  REQUIRE(&elem0.vertex(6) == &elem1.vertex(6));
  REQUIRE(&elem0.vertex(2) != &elem1.vertex(6)); // again, basic sanity check that they're not just all the same
}

TEST_CASE("Refined_connection<Deformed_element>")
{
  cartdg::Storage_params params {3, 5, 3, 6};
  cartdg::Deformed_element coarse (params);
  cartdg::Deformed_element elem0 (params);
  cartdg::Deformed_element elem1 (params);
  cartdg::Deformed_element elem2 (params);
  cartdg::Deformed_element elem3 (params);
  std::vector<cartdg::Deformed_element*> elem_ptrs {&elem0, &elem1};
  REQUIRE_THROWS(cartdg::Refined_connection<cartdg::Deformed_element>(&coarse, elem_ptrs, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}));
  elem_ptrs.push_back(&elem2);
  elem_ptrs.push_back(&elem3);
  SECTION("not reversed")
  {
    cartdg::Refined_connection<cartdg::Deformed_element> con {&coarse, elem_ptrs, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}};
    REQUIRE(con.refined_face.coarse_face() == coarse.face() + (2*0 + 1)*5*6*6);
    REQUIRE(con.direction().i_dim[1] == 2);
    REQUIRE(con.direction().face_sign[1] == 1);
    auto& fine_con = con.connection(1);
    REQUIRE(fine_con.direction().i_dim[0] == 0);
    REQUIRE(fine_con.direction().i_dim[1] == 2);
    REQUIRE(&fine_con.element(0) == &coarse);
    REQUIRE(&fine_con.element(1) == &elem2); // note transposed
    REQUIRE(fine_con.face(0) == con.refined_face.fine_face(1));
    REQUIRE(fine_con.face(1) == elem2.face() + (2*2 + 1)*5*6*6);
    REQUIRE(&coarse.vertex(4) == &elem0.vertex(1));
    REQUIRE(&coarse.vertex(5) == &elem2.vertex(5));
    REQUIRE(&coarse.vertex(6) == &elem1.vertex(3));
    REQUIRE(&coarse.vertex(7) == &elem3.vertex(7));
    REQUIRE(&coarse.vertex(5) != &elem0.vertex(5));
    elem0.vertex_time_step_scale(1) = 0.;
    con.matcher.match(&cartdg::Element::vertex_time_step_scale);
    REQUIRE(elem1.vertex_time_step_scale(1) == Approx(0.5));
    REQUIRE(elem2.vertex_time_step_scale(3) == Approx(0.75));
    REQUIRE(elem2.vertex_time_step_scale(2) == Approx(1.));
  }
  SECTION("reversed")
  {
    // the result should be the same as above, but with fine elements as the left side of the connection and coarse as the right
    cartdg::Refined_connection<cartdg::Deformed_element> con {&coarse, elem_ptrs, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true};
    REQUIRE(con.refined_face.coarse_face() == coarse.face() + (2*2 + 1)*5*6*6);
    auto& fine_con = con.connection(1);
    REQUIRE(fine_con.direction().i_dim[0] == 0);
    REQUIRE(fine_con.direction().i_dim[1] == 2);
    REQUIRE(&fine_con.element(0) == &elem1);
    REQUIRE(&fine_con.element(1) == &coarse);
    REQUIRE(fine_con.face(0) == elem1.face() + (2*0 + 1)*5*6*6);
    REQUIRE(fine_con.face(1) == con.refined_face.fine_face(2)); // note transposed
    REQUIRE(&elem0.vertex(4) == &coarse.vertex(1));
    REQUIRE(&elem1.vertex(5) == &coarse.vertex(5));
    REQUIRE(&elem2.vertex(6) == &coarse.vertex(3));
    REQUIRE(&elem3.vertex(7) == &coarse.vertex(7));
    REQUIRE(&elem0.vertex(5) != &coarse.vertex(5));
  }
  SECTION("stretched")
  {
    std::vector<cartdg::Deformed_element*> elems2 {&elem0, &elem1};
    SECTION("invalid construction throws")
    {
      REQUIRE_THROWS(cartdg::Refined_connection<cartdg::Deformed_element>{&coarse, elem_ptrs, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {true, false}});
      REQUIRE_THROWS(cartdg::Refined_connection<cartdg::Deformed_element>{&coarse, elem_ptrs, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {false, true}});
      REQUIRE_THROWS(cartdg::Refined_connection<cartdg::Deformed_element>{&coarse, elem_ptrs, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {true, true}});
      REQUIRE_THROWS(cartdg::Refined_connection<cartdg::Deformed_element>{&coarse, elems2, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {false, false}});
      REQUIRE_THROWS(cartdg::Refined_connection<cartdg::Deformed_element>{&coarse, elems2, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {true, true}});
    }
    SECTION("stretch dimension 0")
    {
      cartdg::Refined_connection<cartdg::Deformed_element> con {&coarse, elems2, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {true, false}};
      REQUIRE(con.refined_face.coarse_face() == coarse.face() + (2*2 + 1)*5*6*6);
      auto& fine_con = con.connection(1);
      REQUIRE(fine_con.direction().i_dim[0] == 0);
      REQUIRE(fine_con.direction().i_dim[1] == 2);
      REQUIRE(&fine_con.element(0) == &elem1);
      REQUIRE(&fine_con.element(1) == &coarse);
      REQUIRE(fine_con.face(0) == elem1.face() + (2*0 + 1)*5*6*6);
      REQUIRE(fine_con.face(1) == con.refined_face.fine_face(1));
      REQUIRE(&elem0.vertex(4) == &coarse.vertex(1));
      REQUIRE(&elem1.vertex(5) == &coarse.vertex(5));
      REQUIRE(&elem0.vertex(6) == &coarse.vertex(3));
      REQUIRE(&elem1.vertex(7) == &coarse.vertex(7));
      REQUIRE(&elem0.vertex(5) != &coarse.vertex(5));
    }
    SECTION("stretch dimension 1")
    {
      cartdg::Refined_connection<cartdg::Deformed_element> con {&coarse, elems2, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {false, true}};
      REQUIRE(con.refined_face.coarse_face() == coarse.face() + (2*2 + 1)*5*6*6);
      auto& fine_con = con.connection(1);
      REQUIRE(fine_con.direction().i_dim[0] == 0);
      REQUIRE(fine_con.direction().i_dim[1] == 2);
      REQUIRE(&fine_con.element(0) == &elem1);
      REQUIRE(&fine_con.element(1) == &coarse);
      REQUIRE(fine_con.face(0) == elem1.face() + (2*0 + 1)*5*6*6);
      REQUIRE(fine_con.face(1) == con.refined_face.fine_face(1));
      REQUIRE(&elem0.vertex(4) == &coarse.vertex(1));
      REQUIRE(&elem0.vertex(5) == &coarse.vertex(5));
      REQUIRE(&elem1.vertex(6) == &coarse.vertex(3));
      REQUIRE(&elem1.vertex(7) == &coarse.vertex(7));
      REQUIRE(&elem1.vertex(5) != &coarse.vertex(5));
    }
    SECTION("stretch both")
    {
      cartdg::Refined_connection<cartdg::Deformed_element> con {&coarse, {&elem0}, cartdg::Con_dir<cartdg::Deformed_element>{{0, 2}, {1, 1}}, true, {true, true}};
      REQUIRE(con.refined_face.coarse_face() == coarse.face() + (2*2 + 1)*5*6*6);
      auto& fine_con = con.connection(0);
      REQUIRE(fine_con.direction().i_dim[0] == 0);
      REQUIRE(fine_con.direction().i_dim[1] == 2);
      REQUIRE(&fine_con.element(0) == &elem0);
      REQUIRE(&fine_con.element(1) == &coarse);
      REQUIRE(fine_con.face(0) == elem0.face() + (2*0 + 1)*5*6*6);
      REQUIRE(fine_con.face(1) == con.refined_face.fine_face(0));
      REQUIRE(&elem0.vertex(4) == &coarse.vertex(1));
      REQUIRE(&elem0.vertex(5) == &coarse.vertex(5));
      REQUIRE(&elem0.vertex(6) == &coarse.vertex(3));
      REQUIRE(&elem0.vertex(7) == &coarse.vertex(7));
      REQUIRE(&elem0.vertex(5) != &coarse.vertex(1));
    }
  }
}
