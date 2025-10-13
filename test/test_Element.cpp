#include <catch2/catch_all.hpp>
#include <hexed/Element.hpp>
#include <hexed/Equidistant.hpp>
#include <hexed/Tree.hpp>
#include <hexed/Gauss_lobatto.hpp>
#include "testing_utils.hpp"

TEST_CASE("Element") {
  hexed::Storage_params params {4, 5, 3, 6};
  int n_dof = params.n_dof();
  hexed::Tree tree(3, 1.);
  hexed::Element element {params, tree};
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
  for (int i_dim = 0; i_dim < 3; ++i_dim) REQUIRE(element.desired_refinement(i_dim) == 0);

  SECTION("is_sharp") {
    hexed::Gauss_lobatto basis(params.row_size);
    SECTION("3d") {
      hexed::next::Mesh_blocks blocks(3, basis);
      element.create_shape(blocks, 3);
      element.create_fake(blocks);
      // coordinates of this tree are irrelevant, we just need a grafted tree so elem1.is_extruded() will be true
      hexed::Tree& tree1 = *tree.graft(hexed::Array<int>::make(0, 0, 0), hexed::Array<hexed::Int>::make(0, 1, 0));
      hexed::Element elem1(params, tree1);
      elem1.create_shape(blocks);
      elem1.create_fake(blocks);
      elem1.glue_shape(element, {std::vector<double>{0., .5, .5}, std::vector<double>{.5, 1., 1.}});
      REQUIRE(elem1.is_extruded());
      REQUIRE(elem1.has_wall());
      SECTION("snapped edge") {
        element.active_shape().boundary_face_3d()->edge(3).snapped_edge = 0;
        REQUIRE(elem1.is_sharp(0) == false);
        REQUIRE(elem1.is_sharp(1) == false);
        REQUIRE(elem1.is_sharp(2) == true);
        element.active_shape().boundary_face_3d()->edge(1).snapped_edge = 2;
        REQUIRE(elem1.is_sharp(0) == false);
        REQUIRE(elem1.is_sharp(1) == false);
        REQUIRE(elem1.is_sharp(2) == true);
        element.active_shape().boundary_face_3d()->edge(0).snapped_edge = 2;
        REQUIRE(elem1.is_sharp(0) == true);
        REQUIRE(elem1.is_sharp(1) == false);
        REQUIRE(elem1.is_sharp(2) == true);
      }
      SECTION("snapped endpoint") {
        element.active_shape().boundary_block()->vertices()[2]->snapped_endpoint = 1;
        REQUIRE(elem1.is_sharp(0) == false);
        REQUIRE(elem1.is_sharp(1) == false);
        REQUIRE(elem1.is_sharp(2) == false);
        element.active_shape().boundary_block()->vertices()[1]->snapped_endpoint = 1;
        REQUIRE(elem1.is_sharp(0) == true);
        REQUIRE(elem1.is_sharp(1) == false);
        REQUIRE(elem1.is_sharp(2) == true);
      }
    }
    SECTION("2d") {
      hexed::Storage_params par2(1, 4, 2, 6);
      hexed::next::Mesh_blocks blocks(2, basis);
      hexed::Tree tree_2d(2, 1.);
      hexed::Element elem(par2, tree);
      hexed::Tree& tree1 = *tree_2d.graft(hexed::Array<int>::make(0, 0), hexed::Array<hexed::Int>::make(1, 0));
      hexed::Element elem1(par2, tree1);
      elem.create_shape(blocks, 1);
      elem.create_fake(blocks);
      elem1.create_shape(blocks);
      elem1.create_fake(blocks);
      elem1.glue_shape(elem, {std::vector<double>{.5, .5}, std::vector<double>{1., 1.}});
      REQUIRE(elem1.is_extruded());
      REQUIRE(elem1.has_wall());
      elem.active_shape().vertex(2).snapped_point = 0;
      REQUIRE(elem1.is_sharp(0) == false);
      REQUIRE(elem1.is_sharp(1) == false);
      elem.active_shape().vertex(3).snapped_point = 0;
      REQUIRE(elem1.is_sharp(0) == false);
      REQUIRE(elem1.is_sharp(1) == true);
    }
  }
}
