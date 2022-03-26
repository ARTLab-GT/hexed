#include <catch2/catch.hpp>
#include <cartdgConfig.hpp>
#include <Accessible_mesh.hpp>

TEST_CASE("Accessible_mesh")
{
  const int row_size = cartdg::config::max_row_size;
  cartdg::Storage_params params {3, 5, 3, row_size};
  cartdg::Accessible_mesh mesh {params, 1.};
  int sn0 = mesh.add_element(0, false, {0, 0});
  int sn1 = mesh.add_element(0, false, {0, 0});
  REQUIRE(sn0 != sn1); // test uniqueness of serial numbers
  int sn2 = mesh.add_element(3,  true, {1, 2});
  REQUIRE(&mesh.element(0, false, sn0) != &mesh.element(0, false, sn1)); // test uniqueness of elements
  REQUIRE_THROWS(&mesh.element(0, false, sn0 + sn1 + 1)); // test that calling on an invalid serial num throws
  REQUIRE_THROWS(&mesh.element(1, false, sn0)); // test that elements are identified with a specific ref level
  REQUIRE_THROWS(&mesh.element(0,  true, sn0)); // test that elements are identified with a specific deformedness
  REQUIRE(mesh.element(3, true, sn2).vertex(0).pos[0] == Approx(1./8.)); // test that ref level and pos are incorporated
  // test sequential access
  REQUIRE(mesh.cartesian().elements().size() == 2);
  REQUIRE(mesh.deformed().elements().size() == 1);
  REQUIRE(mesh.elements().size() == 3);
  // check that each of the elements appears exactly once in `mesh.elements()`
  cartdg::Element* ptrs [] {&mesh.element(0, false, sn0), &mesh.element(0, false, sn1), &mesh.element(3, true, sn2)};
  for (int i_elem = 0; i_elem < 3; ++i_elem) {
    int count = 0;
    for (int j_elem = 0; j_elem < 3; ++j_elem) {
      cartdg::Element& elem {mesh.elements()[j_elem]};
      if (&elem == ptrs[i_elem]) ++count;
    }
    REQUIRE(count == 1);
  }
  // test connections
  int sn3 = mesh.add_element(3, false, {0, 0});
  int sn4 = mesh.add_element(3, true, {1, 2}); // position doesn't really make sense, but I don't really care
  REQUIRE_THROWS(mesh.connect_cartesian(0, {sn0, sn0 + sn1 + 1}, {0})); // connecting non-existent elements throws
  // make some elements for testing refined face connections
  int car0 = mesh.add_element(1, false, {0, 0}); // note: for this purpose, position doesn't matter
  int def0 = mesh.add_element(2, true, {0, 0});
  int def1 = mesh.add_element(2, true, {0, 0});
  int def2 = mesh.add_element(2, true, {0, 0});
  int def3 = mesh.add_element(2, true, {0, 0});

  SECTION("cartesian-cartesian connection")
  {
    mesh.connect_cartesian(0, {sn1, sn0}, {2});
    auto& con = mesh.cartesian().face_connections()[0];
    REQUIRE(con.direction().i_dim == 2);
    REQUIRE(con.face(0) == mesh.element(0, false, sn1).face() + (2*2 + 1)*5*row_size*row_size);
    REQUIRE(con.face(1) == mesh.element(0, false, sn0).face() + (2*2 + 0)*5*row_size*row_size);
    auto& elem_con = mesh.cartesian().element_connections()[0];
    REQUIRE(&elem_con.element(0) == &mesh.element(0, false, sn1));
    REQUIRE(&elem_con.element(1) == &mesh.element(0, false, sn0));
  }

  SECTION("deformed-cartesian connection")
  {
    mesh.connect_cartesian(3, {sn2, sn3}, {1}, {true, false});
    auto& con = mesh.cartesian().face_connections()[0];
    REQUIRE(con.direction().i_dim == 1);
    REQUIRE(con.face(0) == mesh.element(3,  true, sn2).face() + (1*2 + 1)*5*row_size*row_size);
    REQUIRE(con.face(1) == mesh.element(3, false, sn3).face() + (1*2 + 0)*5*row_size*row_size);
  }

  SECTION("refined face connection")
  {
    // check that it can't find elements with the wrong deformedness
    REQUIRE_THROWS(mesh.connect_hanging_cartesian(1, car0, {def0, def1, def2, def3}, {2}, true, false, false));
    // make a good connection
    mesh.connect_hanging_cartesian(1, car0, {def0, def1, def2, def3}, {2}, true, false, true);
    auto& ref_face {mesh.cartesian().refined_faces()[0]};
    REQUIRE(ref_face.coarse_face() == mesh.element(1, false, car0).face() + (2*2 + 1)*5*row_size*row_size);
  }

  SECTION("deformed-deformed connection")
  {
    // if dimension is same, positivity must be different
    REQUIRE_THROWS(mesh.connect_deformed(3, {sn2, sn4}, {{0, 0}, {0, 0}}));
    REQUIRE_THROWS(mesh.connect_deformed(3, {sn2, sn4}, {{0, 0}, {1, 1}}));
    mesh.connect_deformed(3, {sn4, sn2}, {{1, 0}, {0, 1}});
    REQUIRE(mesh.deformed().face_connections().size() == 1);
    auto& con = mesh.deformed().face_connections()[0];
    REQUIRE(con.direction().i_dim[0] == 1);
    REQUIRE(con.direction().i_dim[1] == 0);
    REQUIRE(con.direction().face_sign[0] == 0);
    REQUIRE(con.direction().face_sign[1] == 1);
    REQUIRE(con.face(0) == mesh.element(3, true, sn4).face() + (1*2 + 0)*5*row_size*row_size);
    REQUIRE(con.face(1) == mesh.element(3, true, sn2).face() + (0*2 + 1)*5*row_size*row_size);
  }

  SECTION("boundary conditions")
  {
    int freestream = mesh.add_boundary_condition(new cartdg::Freestream({0, 0, 1., 1e5}));
    int nonpen = mesh.add_boundary_condition(new cartdg::Nonpenetration);
    // check that connecting to an invalid serial number throws
    REQUIRE_THROWS(mesh.connect_boundary(0, 0, sn0, 1, 0, nonpen + freestream + 1));
    SECTION("cartesian")
    {
      mesh.connect_boundary(0, 0, sn1, 1, 0, freestream);
      {
        // check that it got the right face
        auto& bc_cons {mesh.cartesian().boundary_connections()};
        REQUIRE(bc_cons.size() == 1);
        REQUIRE(bc_cons[0].inside_face() == mesh.element(0, 0, sn1).face() + (2*1 + 0)*5*row_size*row_size);
        // check that the boundary connection is in the deformed connection sequence
        auto& cons = mesh.deformed().face_connections();
        REQUIRE(cons.size() == 1);
        REQUIRE(cons[0].face(0) == bc_cons[0].face(0));
        REQUIRE(cons[0].face(1) == bc_cons[0].face(1));
        // ...and not in the Cartesian or element connection sequences
        REQUIRE(mesh.cartesian().face_connections().size() == 0);
        REQUIRE(mesh.element_connections().size() == 0);
      }
      mesh.connect_boundary(0, 0, sn1, 1, 1, freestream);
      mesh.connect_boundary(0, 0, sn0, 0, 1, nonpen);
      {
        auto& bc_cons {mesh.cartesian().boundary_connections()};
        REQUIRE(bc_cons.size() == 3);
        // check that there are two boundary conditions one of which is used twice and the other of
        // which is used once.
        auto& bcs = mesh.boundary_conditions();
        REQUIRE(bcs.size() == 2);
        int uses [2] {};
        for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
          for (int i_bc = 0; i_bc < 2; ++i_bc) {
            if (bc_cons[i_con].boundary_condition() == &bcs[i_bc]) ++uses[i_bc];
          }
        }
        REQUIRE(std::max(uses[0], uses[1]) == 2);
        REQUIRE(std::min(uses[0], uses[1]) == 1);
      }
    }
    SECTION("deformed")
    {
      mesh.connect_boundary(3, true, sn2, 0, 1, nonpen);
      auto& cons {mesh.deformed().face_connections()};
      REQUIRE(cons.size() == 1);
    }
  }

  SECTION("view with multiple connections")
  {
    int sn5 = mesh.add_element(0, false, {0, 0});
    mesh.connect_cartesian(0, {sn1, sn0}, {0});
    mesh.connect_cartesian(0, {sn1, sn5}, {1});
    mesh.connect_hanging_cartesian(1, car0, {def0, def1, def2, def3}, {2}, true, false, true);
    mesh.connect_deformed(3, {sn4, sn2}, {{1, 0}, {0, 1}});
    REQUIRE(mesh.cartesian().face_connections().size() == 6);
    int count0 = 0;
    int count1 = 0;
    auto& cons = mesh.cartesian().face_connections();
    auto& ref_face {mesh.cartesian().refined_faces()[0]};
    for (int i_con = 0; i_con < cons.size(); ++i_con) {
      auto& con = cons[i_con];
      // is this the connection between sn1 and sn0?
      if (   (con.face(0) == mesh.element(0, false, sn1).face() + (0*2 + 1)*5*row_size*row_size)
          && (con.face(1) == mesh.element(0, false, sn0).face() + (0*2 + 0)*5*row_size*row_size)) ++count0;
      // is this the connection between def3 and car0 (or rather, the fine mortar face involved in that connection)?
      if (   (con.face(1) == mesh.element(2, true, def3).face() + (2*2 + 0)*5*row_size*row_size)
          && (con.face(0) == ref_face.fine_face(3))) ++count1;
    }
    REQUIRE(count0 == 1);
    REQUIRE(count1 == 1);
    count0 = 0;
    count1 = 0;
    int count2 = 0;
    auto& elem_cons = mesh.element_connections();
    for (int i_con = 0; i_con < elem_cons.size(); ++i_con) {
      auto& con = elem_cons[i_con];
      if (   (&con.element(0) == &mesh.element(0, false, sn1))
          && (&con.element(1) == &mesh.element(0, false, sn0))) ++count0;
      if (   (&con.element(0) == &mesh.element(1, false, car0))
          && (&con.element(1) == &mesh.element(2,  true, def1))) ++count1;
      if (   (&con.element(0) == &mesh.element(3,  true, sn4))
          && (&con.element(1) == &mesh.element(3,  true, sn2))) ++count2;
    }
    REQUIRE(count0 == 1);
    REQUIRE(count1 == 1);
    REQUIRE(count2 == 1);
    REQUIRE(elem_cons.size() == 7);
  }

  SECTION("connection validity testing")
  {
    cartdg::Storage_params params1 {3, 4, 2, row_size};
    cartdg::Accessible_mesh mesh1 {params1, 0.64};
    // serial numbers:
    int car [4];
    int def [4];
    int* kinds [] {car, def};
    int coarse [2];
    /*
     * create the following grid (elements identified by serial numbers)
     * +------+------+-------------+
     * |      |      |             |
     * |def[1]|def[3]|             |
     * +------+------+             |
     * |      |      |             |
     * |def[0]|def[2]| coarse[1]   |
     * +------+------+-------------+
     * |      |      |             |
     * |car[1]|car[3]|             |
     * +------+------+             |
     * |      |      |             |
     * |car[0]|car[2]| coarse[0]   |
     * +------+------+-------------+
     */
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int kind = 0; kind < 2; ++kind) {
          kinds[kind][2*i + j] = mesh1.add_element(1, kind, {i, j});
        }
      }
      coarse[i] = mesh1.add_element(0, false, {1, i});
    }
    {
      auto con_val = mesh1.valid();
      REQUIRE(con_val.n_duplicate == 0);
      REQUIRE(con_val.n_missing == 4*(4+4+2)); // every face has a missing connection
      REQUIRE(!con_val);
      REQUIRE_THROWS(con_val.assert_valid());
    }
    // connect the two coarse elements
    mesh1.connect_cartesian(0, {coarse[0], coarse[1]}, {1});
    {
      auto con_val = mesh1.valid();
      REQUIRE(con_val.n_missing == 4*(4+4+2) - 2);
    }
    // connect the fine elements
    for (int i = 0; i < 2; ++i) {
      mesh1.connect_cartesian(1, {car[  i], car[  i+2]}, {0});
      mesh1.connect_cartesian(1, {car[2*i], car[2*i+1]}, {1});
      mesh1.connect_cartesian(1, {car[2*i+1], def[2*i]}, {1}, {false, true});
      mesh1.connect_deformed (1, {def[  i], def[  i+2]}, {{0, 0}, {1, 0}});
      mesh1.connect_deformed (1, {def[2*i], def[2*i+1]}, {{1, 1}, {1, 0}});
    }
    for (int kind = 0; kind < 2; ++kind) {
      mesh1.connect_hanging_cartesian(0, coarse[kind], {kinds[kind][2], kinds[kind][3]}, {0}, 0, false, kind);
    }
    // add boundary conditions
    int bcsn = mesh1.add_boundary_condition(new cartdg::Freestream {{0., 0., 1., 1.}});
    for (int i = 0; i < 2; ++i) {
      for (int kind = 0; kind < 2; ++kind) {
        mesh1.connect_boundary(1, kind, kinds[kind][i], 0, 0, bcsn); // left face
        mesh1.connect_boundary(1, kind, kinds[kind][2*i + kind], 1, kind, bcsn); // bottom/top face
      }
      mesh1.connect_boundary(0, false, coarse[i], 0, 1, bcsn); // right face
      mesh1.connect_boundary(0, false, coarse[i], 1, i, bcsn); // bottom/top face
    }
    // everything should be fine now
    {
      auto con_val = mesh1.valid();
      REQUIRE(con_val.n_duplicate == 0);
      REQUIRE(con_val.n_missing == 0); // every face has a missing connection
      REQUIRE(bool(con_val));
      con_val.assert_valid();
    }
    // check that it can detect duplicate connections
    mesh1.connect_cartesian(0, {coarse[0], coarse[1]}, {1});
    mesh1.connect_cartesian(0, {coarse[0], coarse[1]}, {1});
    {
      auto con_val = mesh1.valid();
      REQUIRE(con_val.n_missing == 0);
      REQUIRE(con_val.n_duplicate == 4);
      REQUIRE(!con_val);
      REQUIRE_THROWS(con_val.assert_valid());
    }
  }
  SECTION("vertices")
  {
    // check that the number of vertices is correct
    auto vertices {mesh.vertices()};
    REQUIRE(vertices.size() == 4*mesh.elements().size());
    // spot-check: vertex 2 of element sn1 should be there
    int count = 0;
    for (int i_vert = 0; i_vert < vertices.size(); ++i_vert) {
      if (&vertices[i_vert] == &mesh.element(0, false, sn1).vertex(2)) ++count;
    }
    REQUIRE(count == 1);
    // test that eaten vertices get dropped from the list
    mesh.connect_cartesian(0, {sn0, sn1}, {0});
    REQUIRE(mesh.vertices().size() == 4*mesh.elements().size() - 2);
  }
}
