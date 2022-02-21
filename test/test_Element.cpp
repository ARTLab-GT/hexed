#include <catch2/catch.hpp>
#include <Element.hpp>
#include "testing_utils.hpp"

TEST_CASE("Element")
{
  cartdg::Storage_params params {4, 5, 3, 6};
  int n_dof = params.n_dof();
  cartdg::Element element {params};
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
  for (int i_face_data = 0; i_face_data < n_dof/params.row_size*3*2; ++i_face_data) {
    element.face()[i_face_data] = 1.;
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
  for (int i_face_data = 0; i_face_data < n_dof/params.row_size*3*2; ++i_face_data) {
    REQUIRE(element.face()[i_face_data] == 1.);
  }
  // test that viscosity is initialized to 0
  for (int i_vert = 0; i_vert < 8; ++i_vert) {
    REQUIRE(element.viscosity()[i_vert] == 0.);
  }
  // test that vertex time step scale is initialized to 1
  for (int i_vert = 0; i_vert < 8; ++i_vert) {
    REQUIRE(element.vertex_time_step_scale()[i_vert] == 1.);
  }
  REQUIRE(element.viscous() == false);
  element.viscosity()[3] = 0.1;
  REQUIRE(element.viscous() == true);
  element.derivative()[0] = 1.;
  element.derivative()[6*6*6 - 1] = 1.;
  REQUIRE(element.derivative()[0] == 1.);
  REQUIRE(element.derivative()[6*6*6 - 1] == 1.);
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint)
  {
    REQUIRE(element.jacobian(0, 0, i_qpoint) == 1.);
    REQUIRE(element.jacobian(1, 0, i_qpoint) == 0.);
    REQUIRE(element.jacobian(0, 2, i_qpoint) == 0.);
    REQUIRE(element.jacobian(2, 2, i_qpoint) == 1.);
    REQUIRE(element.jacobian_determinant(i_qpoint) == 1.);
  }

  SECTION("vertex arrangement")
  {
    // check that vertices start out in correct location
    cartdg::Storage_params params3d {3, 5, 3, 4};
    cartdg::Element elem3d {params3d, {1, 2, -1}, 0.2};
    assert_equal(elem3d.vertex(0).pos, {0.2, 0.4, -0.2});
    assert_equal(elem3d.vertex(1).pos, {0.2, 0.4,  0. });
    assert_equal(elem3d.vertex(2).pos, {0.2, 0.6, -0.2});
    assert_equal(elem3d.vertex(5).pos, {0.4, 0.4,  0. });
    assert_equal(elem3d.vertex(7).pos, {0.4, 0.6,  0. });
    REQUIRE( cartdg::Vertex::are_neighbors(elem3d.vertex(0), elem3d.vertex(1)));
    REQUIRE( cartdg::Vertex::are_neighbors(elem3d.vertex(0), elem3d.vertex(2)));
    REQUIRE(!cartdg::Vertex::are_neighbors(elem3d.vertex(0), elem3d.vertex(3)));
  }

  SECTION("push/fetch viscosity")
  {
    // test push_required_visc
    element.viscosity()[0] = 0.1;
    element.viscosity()[1] = 0.;
    element.viscosity()[2] = 0.;
    element.viscosity()[3] = 0.2;
    element.push_shareable_value(&cartdg::Element::viscosity);
    REQUIRE(element.vertex(0).shared_value() == Approx(0.1));
    REQUIRE(element.vertex(1).shared_value() == Approx(0.));
    REQUIRE(element.vertex(3).shared_value() == Approx(0.2));
    // make sure vertex combination doesn't break anything
    cartdg::Vertex::Transferable_ptr ptr ({0, 0, 0});
    ptr->eat(element.vertex(0));
    REQUIRE(element.vertex(0).shared_value() == Approx(0.1));
    ptr.shareable_value = 0.3;
    REQUIRE(element.vertex(0).shared_value() == Approx(0.3));
    // test fetch_visc
    element.fetch_shareable_value(&cartdg::Element::viscosity);
    REQUIRE(element.viscosity()[0] == Approx(0.3));
    REQUIRE(element.viscosity()[1] == Approx(0.));
    REQUIRE(element.viscosity()[3] == Approx(0.2));
  }
}
