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
}
