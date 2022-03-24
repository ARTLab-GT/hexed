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
    REQUIRE(&fine_con.element(1) == &elem1);
    REQUIRE(fine_con.face(0) == con.refined_face.fine_face(1));
    REQUIRE(fine_con.face(1) == elem1.face() + (2*2 + 1)*5*6*6);
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
    REQUIRE(fine_con.face(1) == con.refined_face.fine_face(1));
  }
}

TEST_CASE("vertex_inds")
{
  SECTION("Same direction")
  {
    //grid3.connect({0, 1}, {0, 0}, {1, 0});
    auto inds = cartdg::vertex_inds(3, {{0, 0}, {1, 0}});
    REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 0);
    REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 1);
    REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 2);
    REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 3);

    //grid3.connect({0, 3}, {2, 2}, {0, 1});
    inds = cartdg::vertex_inds(3, {{2, 2}, {0, 1}});
    REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 0);
    REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 1);
    REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 4);
    REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 5);

    //grid3.connect({2, 0}, {1, 1}, {0, 1});
    inds = cartdg::vertex_inds(3, {{1, 1}, {0, 1}});
    REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 1);
    REQUIRE(inds[0][1] == 2); REQUIRE(inds[1][1] == 3);
    REQUIRE(inds[0][2] == 4); REQUIRE(inds[1][2] == 5);
    REQUIRE(inds[0][3] == 6); REQUIRE(inds[1][3] == 7);
  }

  SECTION("Different direction")
  {
    SECTION("0+ 1+")
    {
      //grid3.connect({0, 1}, {0, 1}, {1, 1});
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {1, 1}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }

    SECTION("0- 1-")
    {
      //grid3.connect({0, 1}, {0, 1}, {0, 0});
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {0, 0}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 1);
      REQUIRE(inds[0][2] == 2); REQUIRE(inds[1][2] == 4);
      REQUIRE(inds[0][3] == 3); REQUIRE(inds[1][3] == 5);
    }

    SECTION("2+ 1+")
    {
      //grid3.connect({0, 1}, {2, 1}, {1, 1});
      auto inds = cartdg::vertex_inds(3, {{2, 1}, {1, 1}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }

    SECTION("0+ 2+")
    {
      //grid3.connect({0, 1}, {0, 2}, {1, 1});
      auto inds = cartdg::vertex_inds(3, {{0, 2}, {1, 1}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 1);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 5);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 3);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }

    SECTION("2+ 0+")
    {
      //grid3.connect({1, 0}, {2, 0}, {1, 1});
      auto inds = cartdg::vertex_inds(3, {{2, 0}, {1, 1}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 1);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 5);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 3);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }

    SECTION("0+ 1-")
    {
      //grid3.connect({0, 1}, {0, 1}, {1, 0});
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {1, 0}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 4);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 5);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 1);
    }

    SECTION("0- 1+")
    {
      //grid3.connect({1, 0}, {0, 1}, {0, 1});
      auto inds = cartdg::vertex_inds(3, {{0, 1}, {0, 1}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 1);
    }

    SECTION("1+ 0-")
    {
      //grid3.connect({0, 1}, {1, 0}, {1, 0});
      auto inds = cartdg::vertex_inds(3, {{1, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 1);
    }

    SECTION("1+ 2-")
    {
      //grid3.connect({0, 1}, {1, 2}, {1, 0});
      auto inds = cartdg::vertex_inds(3, {{1, 2}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 0);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 4);
    }

    SECTION("2+ 0-")
    {
      //grid3.connect({0, 1}, {2, 0}, {1, 0});
      auto inds = cartdg::vertex_inds(3, {{2, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 1);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 2);
    }

    SECTION("2D")
    {
      //grid2.connect({0, 1}, {0, 0}, {1, 0});
      auto inds = cartdg::vertex_inds(2, {{0, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 1);

      //grid2.connect({0, 2}, {0, 1}, {0, 1});
      inds = cartdg::vertex_inds(2, {{0, 1}, {0, 1}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 3);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 1);

      //grid2.connect({0, 3}, {1, 0}, {0, 0});
      inds = cartdg::vertex_inds(2, {{1, 0}, {0, 0}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 2); REQUIRE(inds[1][1] == 1);
    }
  }
}
