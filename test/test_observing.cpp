#include <catch.hpp>

#include <get_mcs_cpg_euler.hpp>
#include <Kernel_settings.hpp>

TEST_CASE("Max characteristic speed")
{
  cartdg::Kernel_settings settings;

  SECTION("1D")
  {
    double sound_speed0 = 400;
    double sound_speed1 = 200;
    double mass = 0.9;
    double heat_rat = 1.3;
    double int_ener0 = mass*sound_speed0*sound_speed0/(heat_rat*(heat_rat - 1));
    double int_ener1 = mass*sound_speed1*sound_speed1/(heat_rat*(heat_rat - 1));
    double read [12] {mass*50, mass*10,
                      mass, mass,
                      int_ener1 + 0.5*mass*250, int_ener0 + 0.5*mass*100,
                      mass*-20, 0,
                      mass, mass,
                      int_ener0 + 0.5*mass*400, int_ener1};
    settings.cpg_heat_rat = heat_rat;
    double mcs = cartdg::get_mcs_cpg_euler(1, 2)(read, 2, settings);
    REQUIRE(mcs == Approx(420));
  }
  SECTION("2D")
  {
    double read[4][4];
    for (int i_qpoint = 0; i_qpoint < 4; ++i_qpoint)
    {
      read[0][i_qpoint] = 2.25;
      read[1][i_qpoint] = 24.5;
      read[2][i_qpoint] = 1.225;
      read[3][i_qpoint] = 101235/0.4 + 0.5*1.225*500;
    }
    double mcs = cartdg::get_mcs_cpg_euler(2, 2)(read[0], 1, settings);
    REQUIRE(mcs == Approx(360).epsilon(0.01));
  }
}
