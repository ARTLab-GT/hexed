#include <catch2/catch_all.hpp>
#include <hexed/Element.hpp>
#include <hexed/Equidistant.hpp>
#include "testing_utils.hpp"

TEST_CASE("Element") {
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
    REQUIRE(!element.is_connected(i_face));
  }
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(0)[i_dof] = 1.2;
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) REQUIRE(element.stage(3)[i_dof] == 0.);
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(3)[i_dof] = 1.3;
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) {
    REQUIRE(element.stage(0)[i_dof] == 1.2);
    REQUIRE(element.stage(1)[i_dof] == 0.);
    REQUIRE(element.stage(2)[i_dof] == 0.);
    REQUIRE(element.stage(3)[i_dof] == 1.3);
  }
  // time step scale exists and should be initialized to 1.
  REQUIRE(element.time_step_scale()[0] == 1.);
  REQUIRE(element.time_step_scale()[params.n_qpoint() - 1] == 1.);
  // artificial viscosity coefficient exists and should be initialized to 0.
  REQUIRE(element.bulk_av_coef()[0] == 0.);
  REQUIRE(element.bulk_av_coef()[params.n_qpoint() - 1] == 0.);
  REQUIRE(element.art_visc_forcing()[0] == 0.);
  REQUIRE(element.art_visc_forcing()[4*params.n_qpoint() - 1] == 0.);
  REQUIRE(element.advection_state()[0] == 0.);
  REQUIRE(element.advection_state()[6*params.n_qpoint() - 1] == 0.);
  // test that vertex time step scale is initialized to the nominal cell size divided by the number of dimensions
  for (int i_vert = 0; i_vert < 8; ++i_vert) {
    REQUIRE(element.vertex_time_step_scale(i_vert) == 1./3.);
  }
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
    REQUIRE(element.jacobian(0, 0, i_qpoint) == 1.);
    REQUIRE(element.jacobian(1, 0, i_qpoint) == 0.);
    REQUIRE(element.jacobian(0, 2, i_qpoint) == 0.);
    REQUIRE(element.jacobian(2, 2, i_qpoint) == 1.);
    REQUIRE(element.jacobian_determinant(i_qpoint) == 1.);
  }
}
