#include <catch2/catch_all.hpp>

#include <hexed/config.hpp>
#include <hexed/math.hpp>
#include <hexed/Deformed_element.hpp>
#include <hexed/Equidistant.hpp>
#include <hexed/Gauss_legendre.hpp>
#include "testing_utils.hpp"

TEST_CASE("Deformed_element")
{
  hexed::Storage_params params {2, 2, 2, 4};
  hexed::Deformed_element element {params};

  // test that accessing the data doesn't segfault
  element.reference_level_normals()[0] = 0.;
  element.reference_level_normals()[4*16 - 1] = 0.;
  element.jacobian_determinant()[0] = 0.;
  element.jacobian_determinant()[16 - 1] = 0.;
  double some_variable = 0;
  REQUIRE(element.face_normal(0) == nullptr);
  REQUIRE(element.face_normal(3) == nullptr);
  element.face_normal(0) = &some_variable;
  element.face_normal(3) = &some_variable;

  for (int i_adj = 0; i_adj < 16; ++i_adj) {
    REQUIRE(element.node_adjustments()[i_adj] == 0.);
  }
  element.node_adjustments()[0] = 2.7;
  element.node_adjustments()[4*2*2 - 1] = 0.04;
  REQUIRE(element.node_adjustments()[0] == 2.7);
  REQUIRE(element.node_adjustments()[4*2*2 - 1] == 0.04);
  REQUIRE(element.vertex(1).is_mobile());

  hexed::Storage_params params3d {1, 1, 3, 2};
  hexed::Deformed_element element3d {params3d};

  SECTION("position calculation")
  {
    const int row_size = 3;
    static_assert (row_size <= hexed::config::max_row_size);
    hexed::Equidistant basis {row_size};
    hexed::Storage_params params2 {2, 4, 2, row_size};
    hexed::Deformed_element elem {params2, {0, 0}, 1., 0, hexed::Mat<2>{.03, .02}};
    elem.vertex(3).pos[0] = 0.63;
    elem.node_adjustments()[2*3 + 1] =  0.2;
    elem.node_adjustments()[3*3 + 1] = -0.1;
    REQUIRE(elem.position(basis, 0)[0] == Catch::Approx(0.03));
    REQUIRE(elem.position(basis, 7)[0] == Catch::Approx(0.83));
    REQUIRE(elem.position(basis, 8)[0] == Catch::Approx(0.63));
    REQUIRE(elem.position(basis, 7)[1] == Catch::Approx(0.52));

    REQUIRE(elem.position(basis, 3)[0] == Catch::Approx(0.53 - 0.2*0.2));
    REQUIRE(elem.position(basis, 4)[0] == Catch::Approx(0.43 - 0.2*(0.2 - 0.1)/2));
    REQUIRE(elem.position(basis, 5)[0] == Catch::Approx(0.33 + 0.2*0.1));
    REQUIRE(elem.position(basis, 3)[1] == Catch::Approx(0.02 + 0.2));
    REQUIRE(elem.position(basis, 4)[1] == Catch::Approx(0.52 + (0.2 - 0.1)/2));
    REQUIRE(elem.position(basis, 5)[1] == Catch::Approx(1.02 - 0.1));
    // check that the face quadrature points are the same as the interior quadrature points
    // that happen to lie on the faces (true for equidistant and Lobatto bases but not Legendre)
    REQUIRE(elem.face_position(basis, 0, 2)[1] == elem.position(basis, 2)[1]);
    REQUIRE(elem.face_position(basis, 2, 1)[0] == elem.position(basis, 3)[0]);
    REQUIRE(elem.face_position(basis, 3, 1)[1] == elem.position(basis, 5)[1]);

    hexed::Deformed_element elem1 {params2};
    elem1.node_adjustments()[1] = 0.1;
    REQUIRE(elem1.position(basis, 0)[0] == Catch::Approx(0.0));
    REQUIRE(elem1.position(basis, 6)[0] == Catch::Approx(1.0));
    REQUIRE(elem1.position(basis, 4)[0] == Catch::Approx(0.55));

    hexed::Storage_params params3 {2, 5, 3, row_size};
    hexed::Deformed_element elem2 {params3, {}, 0.2};
    elem2.node_adjustments()[4] = 0.01;
    REQUIRE(elem2.position(basis, 13)[0] == Catch::Approx(0.101));
    REQUIRE(elem2.position(basis, 13)[1] == Catch::Approx(.1));
    REQUIRE(elem2.position(basis, 13)[2] == Catch::Approx(.1));

    hexed::Gauss_legendre leg_basis {row_size};
    hexed::Deformed_element elem3 {params2, {}, 0.2};
    elem3.node_adjustments()[1] = 0.1;
    elem3.node_adjustments()[3] = -0.2;
    REQUIRE(elem3.position(leg_basis, 3)[0] == Catch::Approx(0.08));
    REQUIRE(elem3.position(leg_basis, 4)[0] == Catch::Approx(0.11));
  }

  SECTION("jacobian calculation")
  {
    const int row_size = 3;
    double faces [6][5*row_size*row_size];
    static_assert (row_size <= hexed::config::max_row_size);
    hexed::Equidistant basis {row_size};
    hexed::Storage_params params2 {2, 4, 2, row_size};
    hexed::Deformed_element elem0 {params2, {0, 0}, 0.2};
    hexed::Deformed_element elem1 {params2, {1, 1}, 0.2};
    elem0.vertex(3).pos = {0.8*0.2, 0.8*0.2, 0.};
    elem1.node_adjustments()[6 + 1] = 0.1;
    // jacobian is correct
    for (int i_face = 0; i_face < 6; ++i_face) elem0.set_face(i_face, faces[i_face]);
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
    for (int i_face = 0; i_face < 6; ++i_face) elem1.set_face(i_face, faces[i_face]);
    elem1.set_jacobian(basis);
    REQUIRE(elem1.jacobian(0, 0, 5) == Catch::Approx(1.));
    REQUIRE(elem1.jacobian(0, 1, 5) == Catch::Approx(0.));
    REQUIRE(elem1.jacobian(1, 0, 5) == Catch::Approx(0.));
    REQUIRE(elem1.jacobian(1, 1, 5) == Catch::Approx(0.9));
    // surface normal is written to face data
    elem0.set_jacobian(basis);
    REQUIRE(elem0.face(0, false)[0] == Catch::Approx(1.));
    REQUIRE(elem0.face(0, false)[row_size] == Catch::Approx(0.));
    REQUIRE(elem0.face(3, false)[2] == Catch::Approx(.2));
    REQUIRE(elem0.face(3, false)[row_size + 2] == Catch::Approx(.8));
    // check time step scale
    REQUIRE(elem0.vertex_time_step_scale(0) == .2/2);
    REQUIRE(elem0.vertex_time_step_scale(3) == Catch::Approx(.2/2*(.8*.8 - .2*.2)/std::sqrt(.8*.8 + .2*.2)));

    hexed::Storage_params params3 {2, 5, 3, row_size};
    hexed::Deformed_element elem2 {params3, {0, 0, 0}, 0.2};
    elem2.vertex(7).pos = {0.8*0.2, 0.8*0.2, 0.8*0.2};
    for (int i_face = 0; i_face < 6; ++i_face) elem2.set_face(i_face, faces[i_face]);
    elem2.set_jacobian(basis);
    REQUIRE(elem2.jacobian(0, 0,  0) == 1.);
    REQUIRE(elem2.jacobian(0, 0, 26) == Catch::Approx( 0.8));
    REQUIRE(elem2.jacobian(0, 1, 26) == Catch::Approx(-0.2));
    REQUIRE(elem2.jacobian(0, 2, 26) == Catch::Approx(-0.2));
    REQUIRE(elem2.jacobian(2, 1, 26) == Catch::Approx(-0.2));
    REQUIRE(elem2.jacobian(2, 2, 26) == Catch::Approx( 0.8));
  }
}
