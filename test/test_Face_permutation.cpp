#include <catch2/catch_all.hpp>
#include <hexed/config.hpp>
#include <hexed/Face_permutation.hpp>
#include <hexed/Accessible_mesh.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/kernels.hpp>
#include <hexed/Spatial.hpp>
#include <hexed/pde.hpp>
#include <hexed/vertex_inds.hpp>

void test_mesh(hexed::Accessible_mesh& mesh) {
  // construct a mesh that has every possible connection configuration by creating a single element
  // and then extruding all its faces
  mesh.add_element(0, 1, hexed::Array<hexed::Int>::make_uniform({mesh.storage_params().n_dim}, 0));
  mesh.extrude();
  auto& elems = mesh.elements();
  auto params = elems[0].storage_params();
  const int n_face_qpoint = params.n_qpoint()/params.row_size;
  // set the face data based on the physical position
  // this implies that for every face connection, data for both faces should be equal,
  // so it can be used to check that the ordering is correct
  hexed::Gauss_legendre basis(params.row_size);
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    auto pos = elem.face_position(basis);
    for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
      if (elem.is_connected(i_face)) {
        for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
          for (int i_var = 0; i_var < params.n_var; ++i_var) {
            elem.face(i_face, false)[i_var*n_face_qpoint + i_qpoint]
            = pos(i_face/2)(i_face%2)(i_var%params.n_dim)[i_qpoint];
          }
        }
      }
    }
  }
  // perform face permutation and check that the faces are indeed equal
  auto connections = mesh.neighbor_connections(1);
  const int n_fdof = params.n_dof()/params.row_size;
  for (int i_con = 0; i_con < connections.size(); ++i_con) {
    auto& con = connections[i_con];
    auto fp = hexed::compute_face_permutation(params.n_dim, params.row_size, con.get_direction(),
                                              con.face(1).flow_state()(0).data(), hexed::laminar);
    fp->match_faces();
    for (int i_dof = 0; i_dof < n_fdof; ++i_dof) {
      REQUIRE(con.face(0).flow_state()(0)[i_dof] == Catch::Approx(con.face(1).flow_state()(0)[i_dof]).scale(1.));
    }
    fp->restore();
  }
  // check that the data has been properly restored to its original order
  // by comparing it to the value it was originally set to
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    auto pos = elem.face_position(basis);
    for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
      if (elem.is_connected(i_face)) {
        for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
          for (int i_var = 0; i_var < params.n_var; ++i_var) {
            REQUIRE(elem.face(i_face, false)[i_var*n_face_qpoint + i_qpoint]
                    == Catch::Approx(pos(i_face/2)(i_face%2)(i_var%params.n_dim)[i_qpoint]).scale(1.));
          }
        }
      }
    }
  }
}

TEST_CASE("Face_permutation") {
  SECTION("2d") {
    hexed::Accessible_mesh mesh {{1, 4, 2, hexed::config::max_row_size}, 1., hexed::laminar};
    std::vector<std::shared_ptr<hexed::Flow_bc>> bc_ptrs;
    for (int i = 0; i < 4; ++i) bc_ptrs.push_back(std::make_shared<hexed::Copy>());
    mesh.add_tree(bc_ptrs);
    test_mesh(mesh);
  }
  SECTION("3d") {
    hexed::Accessible_mesh mesh {{1, 5, 3, hexed::config::max_row_size}, 1., hexed::laminar};
    std::vector<std::shared_ptr<hexed::Flow_bc>> bc_ptrs;
    for (int i = 0; i < 6; ++i) bc_ptrs.push_back(std::make_shared<hexed::Copy>());
    mesh.add_tree(bc_ptrs);
    test_mesh(mesh);
  }
  SECTION("rotation") {
    constexpr int rs = hexed::config::max_row_size;
    hexed::Array<double> data({5, rs, rs});
    data = 0;
    for (int i = 0; i < rs; ++i) {
      for (int j = 0; j < rs; ++j) {
        data(0)(i)[j] = i + .1*j;
      }
    }
    SECTION("+1") {
      hexed::Spatial<hexed::pde::Navier_stokes<>::Pde, false>::Face_permutation<3, rs> perm({{0, 0}, {1, 0}, 1}, data.data());
      perm.match_faces();
      for (int i = 0; i < rs; ++i) {
        for (int j = 0; j < rs; ++j) {
          REQUIRE(data(0)(i)[j] == Catch::Approx(.1*(rs - 1 - i) + j).scale(1.));
        }
      }
      perm.restore();
      for (int i = 0; i < rs; ++i) {
        for (int j = 0; j < rs; ++j) {
          REQUIRE(data(0)(i)[j] == Catch::Approx(i + .1*j).scale(1.));
        }
      }
    }
    SECTION("-1") {
      hexed::Spatial<hexed::pde::Navier_stokes<>::Pde, false>::Face_permutation<3, rs> perm({{0, 0}, {1, 0}, -1}, data.data());
      perm.match_faces();
      for (int i = 0; i < rs; ++i) {
        for (int j = 0; j < rs; ++j) {
          REQUIRE(data(0)(i)[j] == Catch::Approx(.1*i + rs - 1 - j).scale(1.));
        }
      }
      perm.restore();
      for (int i = 0; i < rs; ++i) {
        for (int j = 0; j < rs; ++j) {
          REQUIRE(data(0)(i)[j] == Catch::Approx(i + .1*j).scale(1.));
        }
      }
    }
    SECTION("+2") {
      hexed::Spatial<hexed::pde::Navier_stokes<>::Pde, false>::Face_permutation<3, rs> perm({{0, 0}, {1, 0}, 2}, data.data());
      perm.match_faces();
      for (int i = 0; i < rs; ++i) {
        for (int j = 0; j < rs; ++j) {
          REQUIRE(data(0)(i)[j] == Catch::Approx(rs - 1 - i + .1*(rs - 1 - j)).scale(1.));
        }
      }
      perm.restore();
      for (int i = 0; i < rs; ++i) {
        for (int j = 0; j < rs; ++j) {
          REQUIRE(data(0)(i)[j] == Catch::Approx(i + .1*j).scale(1.));
        }
      }
    }
  }
  SECTION("`vertex_inds` compatibility") {
    hexed::Array<double> data({5, 2, 2});
    data = 0;
    for (int i = 0; i < 4; ++i) data(0)[i] = i;
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      for (int j_dim = 0; j_dim < 3; ++j_dim) {
        for (bool i_sign : {0, 1}) {
          for (bool j_sign : {0, 1}) {
            if (i_dim == j_dim && i_sign == j_sign) continue;
            for (int rotate : {-1, 0, 1, 2}) {
              hexed::Connection_direction dir {{i_dim, j_dim}, {i_sign, j_sign}, rotate};
              hexed::Spatial<hexed::pde::Navier_stokes<>::Pde, false>::Face_permutation<3, 2> perm(dir, data.data());
              perm.match_faces();
              auto inds = hexed::face_vertex_inds(3, dir);
              for (int i = 0; i < 4; ++i) {
                REQUIRE(data(0)[i] == Catch::Approx(double(inds[i])).scale(1));
              }
              perm.restore();
            }
          }
        }
      }
    }
  }
}
