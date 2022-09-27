#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Local_av0_deformed.hpp>
#include <hexed/Local_av1_deformed.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Qpoint_func.hpp>
#include <hexed/Write_face.hpp>

TEST_CASE("local deformed artificial viscosity")
{
  const int row_size = hexed::config::max_row_size;
  static_assert(row_size >= 5); // must be able to represent quartic polynomial exactly
  hexed::Storage_params params {2, 5, 3, row_size};
  const int n_qpoint = params.n_qpoint();
  hexed::Gauss_legendre basis(row_size);
  std::vector<std::unique_ptr<hexed::Deformed_element>> elems;
  elems.emplace_back(new hexed::Deformed_element(params, {}, .4, 1));
  auto& elem = *elems.back();
  elem.vertex(1).pos = {.05*.2, -.1*.2, 0.94*.2};
  elem.set_jacobian(basis);
  double face_normal [6][3*n_qpoint];
  for (int i_face = 0; i_face < 6; ++i_face) {
    elem.face_normal(i_face) = face_normal[i_face];
    for (int i_data = 0; i_data < n_qpoint/row_size*3; ++i_data) {
      face_normal[i_face][i_data] = elem.face()[i_face*n_qpoint/row_size*5 + i_data];
    }
  }
  def_elem_view elem_view {elems};

  // test on a rational function. Expect only approximate agreement
  double scale [] {.1, .2, -.5};
  double arg [n_qpoint];
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    arg[i_qpoint] = 2.;
    auto pos = elem.position(basis, i_qpoint);
    for (int i_dim = 0; i_dim < 3; ++i_dim) arg[i_qpoint] += scale[i_dim]*pos[i_dim];
    for (int i_var = 0; i_var < 5; ++i_var) {
      for (int i_stage = 0; i_stage < 2; ++i_stage) {
        elem.stage(i_stage)[i_var*n_qpoint + i_qpoint] = 0.;
      }
    }
    for (int i_stage = 0; i_stage < 2; ++i_stage) {
      elem.stage(i_stage)[3*n_qpoint + i_qpoint] = 1./arg[i_qpoint];
    }
    elem.art_visc_coef()[i_qpoint] = arg[i_qpoint];
  }
  hexed::Vector_view<hexed::Element&, std::unique_ptr<hexed::Deformed_element>, &hexed::ptr_convert<hexed::Element&, std::unique_ptr<hexed::Deformed_element>>> car {elems};
  (hexed::Write_face<3, row_size>(basis))(car); // `Local_av0_deformed` expects faces to contain numerical LDG state
  (hexed::Local_av0_deformed<3, row_size>(basis, .3, 1.))(elem_view);
  (hexed::Local_av1_deformed<3, row_size>(basis, .3, 1.))(elem_view);
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    double correct = 0.;
    for (int i_dim = 0; i_dim < 3; ++i_dim) correct += scale[i_dim]*scale[i_dim];
    correct *= .3/arg[i_qpoint]/arg[i_qpoint];
    REQUIRE(hexed::Physical_update()(elem, basis, i_qpoint, 0)[3] == Approx(correct).scale(1.));
  }
}
