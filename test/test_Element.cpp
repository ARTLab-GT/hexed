#include <catch2/catch_all.hpp>
#include <hexed/Element.hpp>
#include <hexed/Equidistant.hpp>
#include "testing_utils.hpp"

TEST_CASE("Element")
{
  hexed::Storage_params params {4, 5, 3, 6};
  int n_dof = params.n_dof();
  hexed::Element element {params};
  // test that Storage_params are the same
  REQUIRE(element.storage_params().n_stage == params.n_stage);
  REQUIRE(element.storage_params().n_var == params.n_var);
  REQUIRE(element.storage_params().n_dim == params.n_dim);
  REQUIRE(element.storage_params().row_size == params.row_size);
  // sometimes we will write to and read from some data just
  // to be sure the storage is there and doesn't overlap
  for (int i_stage = 0; i_stage < 4; ++i_stage) {
    for (int i_dof = 0; i_dof < n_dof; ++i_dof) {
      element.stage(i_stage)[i_dof] = 0.;
    }
  }
  for (int i_face = 0; i_face < 6; ++i_face) {
    REQUIRE(element.face_record[i_face] == 0);
    REQUIRE(element.faces[i_face] == nullptr);
  }
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(0)[i_dof] = 1.2;
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) REQUIRE(element.stage(3)[i_dof] == 0.);
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(3)[i_dof] = 1.3;
  for (int i_dof = 0; i_dof < n_dof; ++i_dof)
  {
    REQUIRE(element.stage(0)[i_dof] == 1.2);
    REQUIRE(element.stage(1)[i_dof] == 0.);
    REQUIRE(element.stage(2)[i_dof] == 0.);
    REQUIRE(element.stage(3)[i_dof] == 1.3);
  }
  // time step scale exists and should be initialized to 1.
  REQUIRE(element.time_step_scale()[0] == 1.);
  REQUIRE(element.time_step_scale()[params.n_qpoint() - 1] == 1.);
  // artificial viscosity coefficient exists and should be initialized to 0.
  REQUIRE(element.art_visc_coef()[0] == 0.);
  REQUIRE(element.art_visc_coef()[params.n_qpoint() - 1] == 0.);
  REQUIRE(element.art_visc_forcing()[0] == 0.);
  REQUIRE(element.art_visc_forcing()[4*params.n_qpoint() - 1] == 0.);
  REQUIRE(element.advection_state()[0] == 0.);
  REQUIRE(element.advection_state()[6*params.n_qpoint() - 1] == 0.);
  // test that vertex time step scale is initialized to the nominal cell size divided by the number of dimensions
  for (int i_vert = 0; i_vert < 8; ++i_vert) {
    REQUIRE(element.vertex_time_step_scale(i_vert) == 1./3.);
  }
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint)
  {
    REQUIRE(element.jacobian(0, 0, i_qpoint) == 1.);
    REQUIRE(element.jacobian(1, 0, i_qpoint) == 0.);
    REQUIRE(element.jacobian(0, 2, i_qpoint) == 0.);
    REQUIRE(element.jacobian(2, 2, i_qpoint) == 1.);
    REQUIRE(element.jacobian_determinant(i_qpoint) == 1.);
  }
  REQUIRE(!element.vertex(2).is_mobile());

  SECTION("vertex arrangement")
  {
    // check that vertices start out in correct location
    hexed::Storage_params params3d {3, 5, 3, 4};
    hexed::Element elem3d {params3d, {1, 2, -1}, 0.8, 2, hexed::Mat<3>{.01, .02, .03}};
    REQUIRE(elem3d.nominal_size() == Catch::Approx(0.2));
    REQUIRE_THAT(elem3d.vertex(0).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{0.21, 0.42, -0.17}, hexed::math::Approx_equal(0, 1e-12)));
    REQUIRE_THAT(elem3d.vertex(1).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{0.21, 0.42,  0.03}, hexed::math::Approx_equal(0, 1e-12)));
    REQUIRE_THAT(elem3d.vertex(2).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{0.21, 0.62, -0.17}, hexed::math::Approx_equal(0, 1e-12)));
    REQUIRE_THAT(elem3d.vertex(5).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{0.41, 0.42,  0.03}, hexed::math::Approx_equal(0, 1e-12)));
    REQUIRE_THAT(elem3d.vertex(7).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{0.41, 0.62,  0.03}, hexed::math::Approx_equal(0, 1e-12)));
    REQUIRE( hexed::Vertex::are_neighbors(elem3d.vertex(0), elem3d.vertex(1)));
    REQUIRE( hexed::Vertex::are_neighbors(elem3d.vertex(0), elem3d.vertex(2)));
    REQUIRE(!hexed::Vertex::are_neighbors(elem3d.vertex(0), elem3d.vertex(3)));
  }

  SECTION("push/fetch viscosity")
  {
    // test push_required_visc
    element.vertex_time_step_scale(0) = 0.1;
    element.vertex_time_step_scale(1) = 0.;
    element.vertex_time_step_scale(2) = 0.;
    element.vertex_time_step_scale(3) = 0.2;
    element.push_shareable_value(&hexed::Element::vertex_time_step_scale);
    REQUIRE(element.vertex(0).shared_value() == Catch::Approx(0.1));
    REQUIRE(element.vertex(1).shared_value() == Catch::Approx(0.));
    REQUIRE(element.vertex(3).shared_value() == Catch::Approx(0.2));
    // make sure vertex combination doesn't break anything
    hexed::Vertex::Transferable_ptr ptr ({0, 0, 0});
    ptr->eat(element.vertex(0));
    REQUIRE(element.vertex(0).shared_value() == Catch::Approx(0.1));
    ptr.shareable_value = 0.3;
    REQUIRE(element.vertex(0).shared_value() == Catch::Approx(0.3));
    // test fetch_visc
    element.fetch_shareable_value(&hexed::Element::vertex_time_step_scale);
    REQUIRE(element.vertex_time_step_scale(0) == Catch::Approx(0.3));
    REQUIRE(element.vertex_time_step_scale(1) == Catch::Approx(0.));
    REQUIRE(element.vertex_time_step_scale(3) == Catch::Approx(0.2));
  }

  SECTION("position calculation")
  {
    const int row_size = 5;
    hexed::Storage_params params {2, 4, 2, row_size};
    hexed::Element elem {params, {1, 2}, 0.31, 0, hexed::Mat<2>{.001, .002}};
    hexed::Equidistant basis {row_size};
    // test first and last qpoints
    REQUIRE(elem.position(basis, 0).size() == 2);
    REQUIRE(elem.position(basis, 0)[0] == Catch::Approx(1*0.31 + .001).scale(1.));
    REQUIRE(elem.position(basis, 0)[1] == Catch::Approx(2*0.31 + .002).scale(1.));
    REQUIRE(elem.position(basis, params.n_qpoint() - 1)[0] == Catch::Approx(2*0.31 + .001).scale(1.));
    REQUIRE(elem.position(basis, params.n_qpoint() - 1)[1] == Catch::Approx(3*0.31 + .002).scale(1.));
    static_assert (row_size%2 == 1); // `row_size` must be odd for the following tests to work
    // test the qpoint at the midpoint of the positive-dimension0 face (the right-hand face)
    REQUIRE(elem.position(basis, row_size*(row_size - 1) + row_size/2)[0] == Catch::Approx(  2*0.31 + .001).scale(1.));
    REQUIRE(elem.position(basis, row_size*(row_size - 1) + row_size/2)[1] == Catch::Approx(2.5*0.31 + .002).scale(1.));
    // test the qpoint at the middle of the element (the mean of the vertex positions)
    REQUIRE(elem.position(basis, params.n_qpoint()/2)[0] == Catch::Approx(1.5*0.31 + .001).scale(1.));
    REQUIRE(elem.position(basis, params.n_qpoint()/2)[1] == Catch::Approx(2.5*0.31 + .002).scale(1.));
    // test face position
    REQUIRE(elem.face_position(basis, 1, 3)[0] == Catch::Approx(0.31*(1. + 1.)   + .001));
    REQUIRE(elem.face_position(basis, 1, 3)[1] == Catch::Approx(0.31*(0.75 + 2.) + .002));

    // make sure that face position works in 3d
    hexed::Storage_params params3 {2, 5, 3, row_size};
    hexed::Element elem3 {params3, {0, 0}, 1., 0, hexed::Mat<3>{.003, .003, .003}};
    REQUIRE(elem3.face_position(basis, 0, 7)[0] == Catch::Approx(0.003));
    REQUIRE(elem3.face_position(basis, 0, 7)[1] == Catch::Approx(0.253));
    REQUIRE(elem3.face_position(basis, 0, 7)[2] == Catch::Approx(0.503));
    REQUIRE(elem3.face_position(basis, 3, 7)[0] == Catch::Approx(0.253));
    REQUIRE(elem3.face_position(basis, 3, 7)[1] == Catch::Approx(1.003));
    REQUIRE(elem3.face_position(basis, 3, 7)[2] == Catch::Approx(0.503));
    REQUIRE(elem3.face_position(basis, 4, 7)[0] == Catch::Approx(0.253));
    REQUIRE(elem3.face_position(basis, 4, 7)[1] == Catch::Approx(0.503));
    REQUIRE(elem3.face_position(basis, 4, 7)[2] == Catch::Approx(0.003));
  }

  SECTION("set_jacobian")
  {
    hexed::Equidistant basis(6);
    double faces[6][5*36];
    for (int i_face = 0; i_face < 6; ++i_face) element.faces[i_face] = faces[i_face];
    element.set_jacobian(basis);
    REQUIRE(element.faces[0][0] == Catch::Approx(1.));
    REQUIRE(element.faces[0][1] == Catch::Approx(1.));
    REQUIRE(element.faces[0][36] == Catch::Approx(0.));
    REQUIRE(element.faces[1][0] == Catch::Approx(1.));
    REQUIRE(element.faces[2][0] == Catch::Approx(0.));
    REQUIRE(element.faces[2][1*36] == Catch::Approx(1.));
    REQUIRE(element.faces[2][2*36] == Catch::Approx(0.));

    // test axisymmetric radius
    hexed::Storage_params params {2, 4, 2, 6, true};
    hexed::Element elem {params, {1, 2}, 0.31, 0, hexed::Mat<2>{.001, .002}};
    for (int i_face = 0; i_face < 4; ++i_face) elem.faces[i_face] = faces[i_face];
    REQUIRE(elem.storage_params().axisymmetric == true);
    elem.set_jacobian(basis);
    REQUIRE(elem.radius()[params.n_qpoint() - 1] == elem.position(basis, params.n_qpoint() - 1)[1]);
  }
}
