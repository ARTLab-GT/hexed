#include <catch2/catch_all.hpp>

#include <hexed/config.hpp>
#include <hexed/math.hpp>
#include <hexed/Deformed_element.hpp>
#include <hexed/Equidistant.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Gauss_lobatto.hpp>
#include "testing_utils.hpp"

TEST_CASE("Deformed_element") {
  hexed::Storage_params params {2, 2, 2, 4};
  hexed::Deformed_element element {params};

  // test that accessing the data doesn't segfault
  element.reference_level_normals()[0] = 0.;
  element.reference_level_normals()[4*16 - 1] = 0.;
  element.jacobian_determinant()[0] = 0.;
  element.jacobian_determinant()[16 - 1] = 0.;

  const int row_size = 3;
  static_assert (row_size <= hexed::config::max_row_size);
  hexed::Gauss_lobatto block_basis(row_size);
  hexed::Equidistant basis {row_size};
  hexed::next::Mesh_blocks blocks2d(2, block_basis);
  hexed::next::Mesh_blocks blocks3d(3, block_basis);
  hexed::Storage_params params2 {2, 4, 2, row_size};
  hexed::Storage_params params3 {2, 5, 3, row_size};

  SECTION("position calculation") {
    hexed::Deformed_element elem {params2, {0, 0}, 1., 0, hexed::Mat<2>{.03, .02}};
    elem.create_shape(blocks2d);
    auto p = elem.shape().vertex(3).point({});
    p[0] = 0.63;
    elem.shape().vertex(3).set_pos(p);
    auto pos = elem.position(basis);
    REQUIRE(pos(0)[0] == Catch::Approx(0.03));
    REQUIRE(pos(0)[7] == Catch::Approx(0.83));
    REQUIRE(pos(0)[8] == Catch::Approx(0.63));
    REQUIRE(pos(1)[7] == Catch::Approx(0.52));

    // check that the face quadrature points are the same as the interior quadrature points
    // that happen to lie on the faces (true for equidistant and Lobatto bases but not Legendre)
    auto face_pos = elem.face_position(basis);
    REQUIRE(face_pos(0)(0)(1)[2] == pos(1)[2]);
    REQUIRE(face_pos(1)(0)(0)[1] == pos(0)[3]);
    REQUIRE(face_pos(1)(1)(1)[1] == pos(1)[5]);

    hexed::Gauss_legendre leg_basis {row_size};
    hexed::Deformed_element elem1 {params3, {}, 0.2};
    SECTION("dimensionality must match") {
      REQUIRE_THROWS(elem1.create_shape(blocks2d));
    }
    elem1.create_shape(blocks3d);
    auto pos1 = elem1.position(leg_basis);
    REQUIRE(pos1(2)[0] == Catch::Approx(.2*leg_basis.node(0)));
    REQUIRE(pos1(2)[2] == Catch::Approx(.2*leg_basis.node(2)));
    REQUIRE(pos1(2)[3] == Catch::Approx(.2*leg_basis.node(0)));
    REQUIRE(pos1(1)[3] == Catch::Approx(.1));
  }

  SECTION("splitting") {
    hexed::Deformed_element elem0(params2, {0, 0}, 1., 0, hexed::Mat<2>{.01, .02});
    hexed::Deformed_element elem1(params2, {0, 0}, 1., 0, hexed::Mat<2>{.01, .02});
    elem0.create_shape(blocks2d, hexed::next::Mesh_blocks::no_face);
    elem1.create_shape(blocks2d, hexed::next::Mesh_blocks::no_face);
    elem0.create_fake(blocks2d);
    elem1.split_shape(elem0, .1, 3);
    auto pos0 {elem0.position(basis)};
    auto pos1 {elem1.position(basis)};
    REQUIRE(pos0(0)[0] == Catch::Approx(0.01));
    REQUIRE(pos0(0)[8] == Catch::Approx(1.01));
    REQUIRE(pos0(1)[0] == Catch::Approx(0.02));
    REQUIRE(pos0(1)[8] == Catch::Approx(0.92));
    REQUIRE(pos1(0)[0] == Catch::Approx(0.01));
    REQUIRE(pos1(0)[8] == Catch::Approx(1.01));
    REQUIRE(pos1(1)[0] == Catch::Approx(0.92));
    REQUIRE(pos1(1)[8] == Catch::Approx(1.02));
  }

  SECTION("jacobian calculation") {
    hexed::Deformed_element elem0 {params2, {0, 0}, 0.2};
    hexed::Deformed_element elem1 {params2, {1, 1}, 0.2};
    elem0.create_shape(blocks2d);
    elem0.shape().vertex(3).set_pos(hexed::Mat<3>{0.8*0.2, 0.8*0.2, 0.});
    elem1.create_shape(blocks2d, 2);
    blocks2d.edges_2d()[0].interior()(0)[1] += .1*.2;
    // jacobian is correct
    elem0.set_jacobian(basis);
    REQUIRE(elem0.jacobian(0, 0, 0) == Catch::Approx(1.));
    REQUIRE(elem0.jacobian(0, 1, 0) == Catch::Approx(0.));
    REQUIRE(elem0.jacobian(1, 0, 0) == Catch::Approx(0.));
    REQUIRE(elem0.jacobian(1, 1, 0) == Catch::Approx(1.));
    REQUIRE(elem0.jacobian(0, 0, 6) == Catch::Approx(1.));
    REQUIRE(elem0.jacobian(0, 1, 6) == Catch::Approx(-0.2));
    REQUIRE(elem0.jacobian(1, 0, 6) == Catch::Approx(0.));
    REQUIRE(elem0.jacobian(1, 1, 6) == Catch::Approx(0.8));
    REQUIRE(elem0.jacobian(0, 0, 8) == Catch::Approx(0.8));
    REQUIRE(elem0.jacobian(0, 1, 8) == Catch::Approx(-0.2));
    REQUIRE(elem0.jacobian(1, 0, 8) == Catch::Approx(-0.2));
    REQUIRE(elem0.jacobian(1, 1, 8) == Catch::Approx(0.8));
    REQUIRE(elem0.jacobian_determinant(6) == Catch::Approx(.8));
    elem1.set_jacobian(basis);
    REQUIRE(elem1.jacobian(0, 0, 5) == Catch::Approx(1.));
    REQUIRE(elem1.jacobian(0, 1, 5) == Catch::Approx(0.));
    REQUIRE(elem1.jacobian(1, 0, 5) == Catch::Approx(0.));
    REQUIRE(elem1.jacobian(1, 1, 5) == Catch::Approx(0.9));
    // surface normal is written to face data
    elem0.set_jacobian(basis);
    REQUIRE(elem0.face(0).normal()[0] == Catch::Approx(1.));
    REQUIRE(elem0.face(0).normal()[row_size] == Catch::Approx(0.));
    REQUIRE(elem0.face(3).normal()[2] == Catch::Approx(.2));
    REQUIRE(elem0.face(3).normal()[row_size + 2] == Catch::Approx(.8));
    // check time step scale
    REQUIRE(elem0.vertex_time_step_scale(0) == .2/2);
    REQUIRE(elem0.vertex_time_step_scale(3) == Catch::Approx(.2/2*(.8*.8 - .2*.2)/std::sqrt(.8*.8 + .2*.2)));

    hexed::Deformed_element elem2 {params3, {0, 0, 0}, 0.2};
    elem2.create_shape(blocks3d);
    elem2.shape().vertex(7).set_pos(hexed::Mat<3>{0.8*0.2, 0.8*0.2, 0.8*0.2});
    elem2.set_jacobian(basis);
    REQUIRE(elem2.jacobian(0, 0,  0) == 1.);
    REQUIRE(elem2.jacobian(0, 0, 26) == Catch::Approx( 0.8));
    REQUIRE(elem2.jacobian(0, 1, 26) == Catch::Approx(-0.2));
    REQUIRE(elem2.jacobian(0, 2, 26) == Catch::Approx(-0.2));
    REQUIRE(elem2.jacobian(2, 1, 26) == Catch::Approx(-0.2));
    REQUIRE(elem2.jacobian(2, 2, 26) == Catch::Approx( 0.8));
  }
}
