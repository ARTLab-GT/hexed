#include <catch2/catch_all.hpp>
#include <hexed/connection.hpp>

TEST_CASE("direction conversion") {
  hexed::Con_dir<hexed::Element> car {2};
  hexed::Con_dir<hexed::Deformed_element> def = car;
  REQUIRE(def.i_dim[0] == 2);
  REQUIRE(def.i_dim[1] == 2);
  REQUIRE(def.face_sign[0] == 1);
  REQUIRE(def.face_sign[1] == 0);
}

TEST_CASE("Element_face_connection<Element>") {
  hexed::Storage_params params {3, 5, 3, 6};
  hexed::Element elem0 (params);
  hexed::Element elem1 (params);
  REQUIRE(!elem0.is_connected(3));
  REQUIRE(!elem1.is_connected(2));
  {
    hexed::Element_face_connection<hexed::Element> con ({&elem0, &elem1}, hexed::Con_dir<hexed::Element>{1});
    REQUIRE(elem0.face(3, false) == con.state(0, false));
    REQUIRE(elem1.face(2, false) == con.state(1, false));
    REQUIRE(con.direction().i_dim == 1);
    REQUIRE(&con.element(0) == &elem0);
    REQUIRE(&con.element(1) == &elem1);
    REQUIRE(con.state(0, false) == elem0.face(3, false));
    REQUIRE(con.state(1, false) == elem1.face(2, false));
  }
  // check that faces get reset to nullptr after the connection is deleted
  REQUIRE(!elem0.is_connected(3));
  REQUIRE(!elem1.is_connected(2));
}

TEST_CASE("Element_face_connection<Deformed_element>") {
  hexed::Storage_params params {3, 5, 3, 6};
  hexed::Deformed_element elem0 (params);
  hexed::Deformed_element elem1 (params);
  hexed::Element_face_connection<hexed::Deformed_element> con ({&elem0, &elem1}, hexed::Con_dir<hexed::Deformed_element>{{2, 1}, {0, 1}});
  REQUIRE(elem0.face(4, false) == con.state(0, false));
  REQUIRE(elem1.face(3, false) == con.state(1, false));
  REQUIRE(con.direction().i_dim[0] == 2);
  REQUIRE(con.direction().i_dim[1] == 1);
  REQUIRE(con.direction().face_sign[0] == 0);
  REQUIRE(con.direction().face_sign[1] == 1);
  REQUIRE(&con.element(0) == &elem0);
  REQUIRE(&con.element(1) == &elem1);
  con.normal(0)[3*6*6 - 1] = 1; // check that normal storage is big enough (otherwise segfault)
  con.normal(1)[3*6*6 - 1] = 1;
  REQUIRE(con.normal() == con.normal(0));
  REQUIRE(con.state(0, false) == elem0.face(4, false));
  REQUIRE(con.state(1, false) == elem1.face(3, false));
}

TEST_CASE("Refined_connection<Deformed_element>") {
  hexed::Storage_params params {3, 5, 3, 6};
  hexed::Deformed_element coarse (params);
  hexed::Deformed_element elem0 (params);
  hexed::Deformed_element elem1 (params);
  hexed::Deformed_element elem2 (params);
  hexed::Deformed_element elem3 (params);
  std::vector<hexed::Deformed_element*> elem_ptrs {&elem0, &elem1};
  REQUIRE_THROWS(hexed::Refined_connection<hexed::Deformed_element>(&coarse, elem_ptrs, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}));
  elem_ptrs.push_back(&elem2);
  elem_ptrs.push_back(&elem3);
  SECTION("not reversed") {
    {
      hexed::Refined_connection<hexed::Deformed_element> con {&coarse, elem_ptrs, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}};
      REQUIRE(con.refined_face.coarse == coarse.face(2*0 + 1, false));
      REQUIRE(con.direction().i_dim[1] == 2);
      REQUIRE(con.direction().face_sign[1] == 1);
      auto& fine_con = con.connection(1);
      REQUIRE(fine_con.direction().i_dim[0] == 0);
      REQUIRE(fine_con.direction().i_dim[1] == 2);
      REQUIRE(&fine_con.element(0) == &coarse);
      REQUIRE(&fine_con.element(1) == &elem2); // note transposed
      REQUIRE(fine_con.state(0, false) == con.refined_face.fine[1]);
      REQUIRE(coarse.face(1, false) == con.coarse_state());
      for (int i_con = 0; i_con < 4; ++i_con) {
        auto& c = con.connection(i_con);
        REQUIRE(c.element(1).face(5, false) == c.state(1, false));
      }
    }
    REQUIRE(!coarse.is_connected(1));
    for (auto ptr : elem_ptrs) {
      REQUIRE(!ptr->is_connected(5));
    }
  }
  SECTION("reversed") {
    // the result should be the same as above, but with fine elements as the left side of the connection and coarse as the right
    hexed::Refined_connection<hexed::Deformed_element> con {&coarse, elem_ptrs, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true};
    REQUIRE(con.refined_face.coarse == coarse.face(2*2 + 1, false));
    auto& fine_con = con.connection(1);
    REQUIRE(fine_con.direction().i_dim[0] == 0);
    REQUIRE(fine_con.direction().i_dim[1] == 2);
    REQUIRE(&fine_con.element(0) == &elem1);
    REQUIRE(&fine_con.element(1) == &coarse);
    REQUIRE(fine_con.state(1, false) == con.refined_face.fine[2]); // note transposed
  }
  SECTION("stretched") {
    std::vector<hexed::Deformed_element*> elems2 {&elem0, &elem1};
    SECTION("invalid construction throws") {
      REQUIRE_THROWS(hexed::Refined_connection<hexed::Deformed_element>{&coarse, elem_ptrs, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true, {true, false}});
      REQUIRE_THROWS(hexed::Refined_connection<hexed::Deformed_element>{&coarse, elem_ptrs, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true, {false, true}});
      REQUIRE_THROWS(hexed::Refined_connection<hexed::Deformed_element>{&coarse, elem_ptrs, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true, {true, true}});
      REQUIRE_THROWS(hexed::Refined_connection<hexed::Deformed_element>{&coarse, elems2, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true, {false, false}});
      REQUIRE_THROWS(hexed::Refined_connection<hexed::Deformed_element>{&coarse, elems2, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true, {true, true}});
    }
    SECTION("stretch dimension 0") {
      hexed::Refined_connection<hexed::Deformed_element> con {&coarse, elems2, hexed::Con_dir<hexed::Deformed_element>{{2, 0}, {1, 1}}, false, {true, false}};
      REQUIRE(con.refined_face.stretch[0] == false);
      REQUIRE(con.refined_face.stretch[1] == true);
      REQUIRE(con.refined_face.coarse == coarse.face(2*2 + 1, false));
      auto& fine_con = con.connection(1);
      REQUIRE(fine_con.direction().i_dim[0] == 2);
      REQUIRE(fine_con.direction().i_dim[1] == 0);
      REQUIRE(&fine_con.element(0) == &coarse);
      REQUIRE(&fine_con.element(1) == &elem1);
      REQUIRE(fine_con.state(0, false) == con.refined_face.fine[1]);
      REQUIRE(fine_con.state(1, false) == elem1.face(2*0 + 1, false));
    }
    SECTION("stretch dimension 0 reverse") {
      hexed::Refined_connection<hexed::Deformed_element> con {&coarse, elems2, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true, {true, false}};
      REQUIRE(con.refined_face.coarse == coarse.face(2*2 + 1, false));
      auto& fine_con = con.connection(1);
      REQUIRE(fine_con.direction().i_dim[0] == 0);
      REQUIRE(fine_con.direction().i_dim[1] == 2);
      REQUIRE(&fine_con.element(0) == &elem1);
      REQUIRE(&fine_con.element(1) == &coarse);
      REQUIRE(fine_con.state(0, false) == elem1.face(2*0 + 1, false));
      REQUIRE(fine_con.state(1, false) == con.refined_face.fine[1]);
    }
    SECTION("stretch dimension 1") {
      hexed::Refined_connection<hexed::Deformed_element> con {&coarse, elems2, hexed::Con_dir<hexed::Deformed_element>{{2, 0}, {1, 1}}, false, {false, true}};
      REQUIRE(con.refined_face.coarse == coarse.face(2*2 + 1, false));
      auto& fine_con = con.connection(1);
      REQUIRE(fine_con.direction().i_dim[0] == 2);
      REQUIRE(fine_con.direction().i_dim[1] == 0);
      REQUIRE(&fine_con.element(0) == &coarse);
      REQUIRE(&fine_con.element(1) == &elem1);
      REQUIRE(fine_con.state(0, false) == con.refined_face.fine[1]);
      REQUIRE(fine_con.state(1, false) == elem1.face(2*0 + 1, false));
    }
    SECTION("stretch both") {
      hexed::Refined_connection<hexed::Deformed_element> con {&coarse, {&elem0}, hexed::Con_dir<hexed::Deformed_element>{{0, 2}, {1, 1}}, true, {true, true}};
      REQUIRE(con.refined_face.coarse == coarse.face(2*2 + 1, false));
      auto& fine_con = con.connection(0);
      REQUIRE(fine_con.direction().i_dim[0] == 0);
      REQUIRE(fine_con.direction().i_dim[1] == 2);
      REQUIRE(&fine_con.element(0) == &elem0);
      REQUIRE(&fine_con.element(1) == &coarse);
      REQUIRE(fine_con.state(0, false) == elem0.face(2*0 + 1, false));
      REQUIRE(fine_con.state(1, false) == con.refined_face.fine[0]);
    }
    SECTION("not transposed") {
      hexed::Refined_connection<hexed::Deformed_element> con {&coarse, {elems2}, hexed::Con_dir<hexed::Deformed_element>{{0, 1}, {1, 1}}, true, {true, false}};
      REQUIRE(con.refined_face.stretch[0] == true);
      REQUIRE(con.refined_face.stretch[1] == false);
    }
  }
}
