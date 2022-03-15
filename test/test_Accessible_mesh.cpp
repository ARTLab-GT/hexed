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
  REQUIRE_THROWS(mesh.connect_cartesian(0, {sn0, sn0 + sn1 + 1}, {0})); // connecting non-existent elements throws
  SECTION("cartesian-cartesian connection")
  {
    mesh.connect_cartesian(0, {sn1, sn0}, {2});
    auto& con = mesh.cartesian().connections()[0];
    REQUIRE(con.direction().i_dim == 2);
    REQUIRE(con.face(0) == mesh.element(0, false, sn1).face() + (2*2 + 1)*5*row_size*row_size);
    REQUIRE(con.face(1) == mesh.element(0, false, sn0).face() + (2*2 + 0)*5*row_size*row_size);
  }
  SECTION("deformed-cartesian connection")
  {
    mesh.connect_cartesian(3, {sn2, sn3}, {1}, {true, false});
    auto& con = mesh.cartesian().connections()[0];
    REQUIRE(con.direction().i_dim == 1);
    REQUIRE(con.face(0) == mesh.element(3,  true, sn2).face() + (1*2 + 1)*5*row_size*row_size);
    REQUIRE(con.face(1) == mesh.element(3, false, sn3).face() + (1*2 + 0)*5*row_size*row_size);
  }
  SECTION("refined face connection")
  {
    int car0 = mesh.add_element(1, false, {0, 0}); // note: for this purpose, position doesn't matter
    int def0 = mesh.add_element(2, true, {0, 0});
    int def1 = mesh.add_element(2, true, {0, 0});
    int def2 = mesh.add_element(2, true, {0, 0});
    int def3 = mesh.add_element(2, true, {0, 0});
    // check that it can't find elements with the wrong deformedness
    REQUIRE_THROWS(mesh.connect_hanging_cartesian(1, car0, {def0, def1, def2, def3}, {2}, true, false, false));
    // make a good connection
    mesh.connect_hanging_cartesian(1, car0, {def0, def1, def2, def3}, {2}, true, false, true);
    auto& ref_face {mesh.cartesian().refined_faces()[0]};
    REQUIRE(ref_face.coarse_face() == mesh.element(1, false, car0).face() + (2*2 + 1)*5*row_size*row_size);
  }
  SECTION("deformed-deformed connection")
  {
    int sn4 = mesh.add_element(3, true, {1, 2}); // position doesn't really make sense, but I don't really care
    // if dimension is same, positivity must be different
    REQUIRE_THROWS(mesh.connect_deformed(3, {sn2, sn4}, {{0, 0}, {0, 0}}));
    REQUIRE_THROWS(mesh.connect_deformed(3, {sn2, sn4}, {{0, 0}, {1, 1}}));
    mesh.connect_deformed(3, {sn4, sn2}, {{1, 0}, {0, 1}});
    auto& con = mesh.deformed().connections()[0];
    REQUIRE(con.direction().i_dim[0] == 1);
    REQUIRE(con.direction().i_dim[1] == 0);
    REQUIRE(con.direction().face_sign[0] == 0);
    REQUIRE(con.direction().face_sign[1] == 1);
    REQUIRE(con.face(0) == mesh.element(3, true, sn4).face() + (1*2 + 0)*5*row_size*row_size);
    REQUIRE(con.face(1) == mesh.element(3, true, sn2).face() + (0*2 + 1)*5*row_size*row_size);
  }
  SECTION("connection vector size")
  {
    mesh.connect_cartesian(0, {sn1, sn0}, {0});
    mesh.connect_cartesian(0, {sn1, sn0}, {1});
    REQUIRE(mesh.cartesian().connections().size() == 2);
  }
}
