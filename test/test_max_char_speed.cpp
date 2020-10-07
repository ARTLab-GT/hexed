#include <catch.hpp>

#include <kernels/max_char_speed/cpg_euler_max.hpp>

TEST_CASE("Max characteristic speed")
{
  SECTION("1D")
  {
    double sound_speed0 = 400;
    double sound_speed1 = 200;
    double mass = 0.9;
    double sp_heat_rat = 1.3;
    double int_ener0 = mass*sound_speed0*sound_speed0/(sp_heat_rat*(sp_heat_rat - 1));
    double int_ener1 = mass*sound_speed1*sound_speed1/(sp_heat_rat*(sp_heat_rat - 1));
    double read [12] {mass*50, mass*10,
                      mass, mass,
                      int_ener1 + 0.5*mass*250, int_ener0 + 0.5*mass*100,
                      mass*-20, 0,
                      mass, mass,
                      int_ener0 + 0.5*mass*400, int_ener1};
    double mcs = cartdg::cpg_euler_max<3, 2, 2>(read, 2, sp_heat_rat);
    REQUIRE(mcs == Approx(420));
  }
  SECTION("2D")
  {
    double read[4] {2.25, 24.5, 1.225, 101235/0.4 + 0.5*1.225*500};
    double mcs = cartdg::cpg_euler_max<4, 1, 1>(read, 1, 1.4);
    REQUIRE(mcs == Approx(360).epsilon(0.01));
  }
}
