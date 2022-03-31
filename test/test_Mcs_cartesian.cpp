#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <Mcs_cartesian.hpp>

TEST_CASE("Mcs_cartesian")
{
  SECTION("1D")
  {
    double sound_speed0 = 400;
    double sound_speed1 = 200;
    double mass = 0.9;
    double heat_rat = 1.3;
    double int_ener0 = mass*sound_speed0*sound_speed0/(heat_rat*(heat_rat - 1));
    double int_ener1 = mass*sound_speed1*sound_speed1/(heat_rat*(heat_rat - 1));
    cartdg::Storage_params params {1, 3, 1, 2};
    cartdg::elem_vec elements;
    elements.emplace_back(new cartdg::Element {params, {}, 1.027});
    elements.emplace_back(new cartdg::Element {params, {}, 1.027});
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
    double mcs = (*cartdg::kernel_factory<cartdg::Mcs_cartesian>(1, 2, heat_rat))(elem_view);
    REQUIRE(mcs == Approx(420*0.6/1.027));
  }
  SECTION("2D")
  {
    cartdg::Storage_params params {3, 4, 2, 2};
    cartdg::elem_vec elements;
    elements.emplace_back(new cartdg::Element {params});
    double* read = elements.back()->stage(0);
    for (int i_qpoint = 0; i_qpoint < 4; ++i_qpoint)
    {
      read[0*4 + i_qpoint] = 2.25;
      read[1*4 + i_qpoint] = 24.5;
      read[2*4 + i_qpoint] = 1.225;
      read[3*4 + i_qpoint] = 101235/0.4 + 0.5*1.225*500;
    }
    car_elem_view elem_view {elements};
    double mcs = (*cartdg::kernel_factory<cartdg::Mcs_cartesian>(2, 2))(elem_view);
    REQUIRE(mcs == Approx(360).epsilon(0.01));
  }
}
