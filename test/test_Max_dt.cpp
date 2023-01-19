#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Spatial.hpp>
#include <hexed/pde.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Max_dt")
{
  hexed::Gauss_legendre basis(2);
  SECTION("cartesian")
  {
    SECTION("1D")
    {
      double sound_speed0 = 400;
      double sound_speed1 = 200;
      double mass = 0.9;
      double heat_rat = 1.4;
      double int_ener0 = mass*sound_speed0*sound_speed0/(heat_rat*(heat_rat - 1));
      double int_ener1 = mass*sound_speed1*sound_speed1/(heat_rat*(heat_rat - 1));
      hexed::Storage_params params {1, 3, 1, 2};
      std::vector<std::unique_ptr<hexed::Element>> elements;
      elements.emplace_back(new hexed::Element {params, {}, 1.027});
      elements.emplace_back(new hexed::Element {params, {}, 1.027});
      double read [2][6] {{mass*50, mass*10,
                           mass, mass,
                           int_ener1 + 0.5*mass*250, int_ener0 + 0.5*mass*100},
                          {mass*-20, 0,
                           mass, mass,
                           int_ener0 + 0.5*mass*400, int_ener1}};
      for (int i_elem : {0, 1}) {
        for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
          elements[i_elem]->time_step_scale()[i_qpoint] = 0.6;
        }
        elements[i_elem]->time_step_scale()[1] = 0.17;
        for (int i_dof = 0; i_dof < 6; ++i_dof) {
          elements[i_elem]->stage(0)[i_dof] = read[i_elem][i_dof];
        }
      }
      car_elem_view elem_view {elements};
      double dt = (*hexed::kernel_factory<hexed::Spatial<hexed::Element, hexed::pde::Navier_stokes<false>::Pde>::Max_dt>(1, 2, basis))(elem_view);
      REQUIRE(dt == Approx(basis.max_cfl_convective()/(420*0.6/1.027)));
    }
    SECTION("2D")
    {
      hexed::Storage_params params {3, 4, 2, 2};
      std::vector<std::unique_ptr<hexed::Element>> elements;
      elements.emplace_back(new hexed::Element {params});
      double* read = elements.back()->stage(0);
      for (int i_qpoint = 0; i_qpoint < 4; ++i_qpoint)
      {
        read[0*4 + i_qpoint] = 2.25;
        read[1*4 + i_qpoint] = 24.5;
        read[2*4 + i_qpoint] = 1.225;
        read[3*4 + i_qpoint] = 101235/0.4 + 0.5*1.225*500;
      }
      car_elem_view elem_view {elements};
      double dt = (*hexed::kernel_factory<hexed::Spatial<hexed::Element, hexed::pde::Navier_stokes<false>::Pde>::Max_dt>(2, 2, basis))(elem_view);
      REQUIRE(dt == Approx(basis.max_cfl_convective()/360./2.).epsilon(0.01));
    }
  }

  SECTION("deformed")
  {
    hexed::Storage_params params {3, 5, 3, 4};
    static_assert (hexed::config::max_row_size >= 4);
    hexed::Gauss_legendre basis(params.row_size);
    int n_qpoint = params.n_qpoint();
    std::vector<std::unique_ptr<hexed::Deformed_element>> elems;
    double faces [6][5*4*4];
    for (int i_elem = 0; i_elem < 3; ++i_elem) {
      elems.emplace_back(new hexed::Deformed_element {params, {}, 0.3});
      double* state = elems.back()->stage(0);
      for (int i_face = 0; i_face < 6; ++i_face) elems.back()->faces[i_face] = faces[i_face];
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
    double dt = (*hexed::kernel_factory<hexed::Spatial<hexed::Deformed_element, hexed::pde::Navier_stokes<false>::Pde>::Max_dt>(3, 4, basis))(elem_view);
    REQUIRE(dt == Approx(basis.max_cfl_convective()/(1e3/0.3)/3.));
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
    dt = (*hexed::kernel_factory<hexed::Spatial<hexed::Deformed_element, hexed::pde::Navier_stokes<false>::Pde>::Max_dt>(3, 4, basis))(elem_view);
    REQUIRE(dt == Approx(basis.max_cfl_convective()/(1e3*0.95*(1. + 2*2.)/3./0.3)/3.));
  }
}
