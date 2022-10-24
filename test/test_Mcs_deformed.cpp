#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Mcs_deformed.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Mcs_deformed")
{
  hexed::Storage_params params {3, 5, 3, 4};
  hexed::Gauss_legendre basis(params.row_size);
  int n_qpoint = params.n_qpoint();
  std::vector<std::unique_ptr<hexed::Deformed_element>> elems;
  for (int i_elem = 0; i_elem < 3; ++i_elem) {
    elems.emplace_back(new hexed::Deformed_element {params, {}, 0.3});
    double* state = elems.back()->stage(0);
    elems.back()->set_jacobian(basis);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      // set state to 0 velocity and speed of sound 100
      for (int i_dim = 0; i_dim < 3; ++i_dim) state[i_dim*n_qpoint + i_qpoint] = 0.;
      state[3*n_qpoint + i_qpoint] = 1.4;
      // set a modestly higher sound speed in element 1 just to mix it up
      state[4*n_qpoint + i_qpoint] = ((i_elem == 1) ? 1e6 : 1e4)/0.4;
    }
  }
  // test that actual characteristic speed is measured correctly
  def_elem_view elem_view {elems};
  REQUIRE((*hexed::kernel_factory<hexed::Mcs_deformed>(3, 4, hexed::char_speed::Inviscid(1.4)))(elem_view) == Approx(1e3/0.3));
  // scale 2 entries of jacobian
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    auto& pos = elems[1]->vertex(i_vert).pos;
    pos[0] *= .5;
    pos[1] *= .5;
  }
  elems[1]->set_jacobian(basis);
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    elems[1]->time_step_scale()[i_qpoint] = 0.9;
  }
  elems[1]->time_step_scale()[2] = 0.95;
  // test that jacobian & time step scale are accounted for
  REQUIRE((*hexed::kernel_factory<hexed::Mcs_deformed>(3, 4, hexed::char_speed::Inviscid(1.4)))(elem_view) == Approx(1e3*0.95*(1. + 2*2.)/3./0.3));
}
