#include <catch2/catch.hpp>
#include <config.hpp>
#include <Face_permutation.hpp>
#include <Accessible_mesh.hpp>
#include <Gauss_legendre.hpp>

void test_mesh(cartdg::Accessible_mesh& mesh)
{
  mesh.add_element(0, 1, {});
  mesh.extrude();
  auto& elems = mesh.elements();
  auto params = elems[0].storage_params();
  const int n_face_qpoint = params.n_qpoint()/params.row_size;
  const int n_face_dof = n_face_qpoint*params.n_var;
  cartdg::Gauss_legendre basis(params.row_size);
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        auto pos = elem.face_position(basis, i_face, i_qpoint);
        for (int i_var = 0; i_var < params.n_var; ++i_var) {
          elem.face()[i_face*n_face_dof + i_var*n_face_qpoint + i_qpoint] = pos[i_var%params.n_dim];
        }
      }
    }
  }
  auto& connections = mesh.deformed().face_connections();
  for (int i_con = 0; i_con < connections.size(); ++i_con) {
    auto& con = connections[i_con];
    auto fp = cartdg::kernel_factory<cartdg::Face_permutation>(params.n_dim, params.row_size, con);
    fp->match_faces();
    for (int i_dof = 0; i_dof < n_face_dof; ++i_dof) {
      REQUIRE(con.face(0)[i_dof] == Approx(con.face(1)[i_dof]).scale(1.));
    }
    fp->restore();
  }
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        auto pos = elem.face_position(basis, i_face, i_qpoint);
        for (int i_var = 0; i_var < params.n_var; ++i_var) {
          REQUIRE(elem.face()[i_face*n_face_dof + i_var*n_face_qpoint + i_qpoint]
                  == Approx(pos[i_var%params.n_dim]).scale(1.));
        }
      }
    }
  }
}

TEST_CASE("Face_permutation")
{
  SECTION("2d") {
    cartdg::Accessible_mesh mesh {{1, 4, 2, cartdg::config::max_row_size}, 1.};
    test_mesh(mesh);
  }
  SECTION("3d") {
    cartdg::Accessible_mesh mesh {{1, 5, 3, cartdg::config::max_row_size}, 1.};
    test_mesh(mesh);
  }
}
