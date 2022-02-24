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
  REQUIRE(mesh.cartesian_elements().size() == 2);
  REQUIRE(mesh.deformed_elements().size() == 1);
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
  REQUIRE_THROWS(mesh.connect_cartesian(0, 0, {sn0, sn0 + sn1 + 1})); // connecting non-existent elements throws
  SECTION("cartesian-cartesian connection")
  {
    mesh.connect_cartesian(0, 2, {sn1, sn0});
    auto con = mesh.cartesian_connections()[0];
    REQUIRE(con.i_dim == 2);
    REQUIRE(con.face[0] == mesh.element(0, false, sn1).face() + (2*2 + 1)*params.n_qpoint());
    REQUIRE(con.face[1] == mesh.element(0, false, sn0).face() + (2*2 + 0)*params.n_qpoint());
  }
  SECTION("deformed-cartesian connection")
  {
    mesh.connect_cartesian(3, 1, {sn2, sn3}, {true, false});
    auto con = mesh.cartesian_connections()[0];
    REQUIRE(con.i_dim == 1);
    REQUIRE(con.face[0] == mesh.element(3,  true, sn2).face() + (1*2 + 1)*params.n_qpoint());
    REQUIRE(con.face[1] == mesh.element(3, false, sn3).face() + (1*2 + 0)*params.n_qpoint());
  }
  SECTION("connection vector size")
  {
    mesh.connect_cartesian(0, 0, {sn1, sn0});
    mesh.connect_cartesian(0, 1, {sn1, sn0});
    REQUIRE(mesh.cartesian_connections().size() == 2);
  }
}
