#include <catch.hpp>

#include <kernels/observing/cpg_euler_max.hpp>
#include <kernels/observing/cpg_euler_physical_step.hpp>

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
    double mcs = cartdg::cpg_euler_max<3, 2, 2>(read, 2, settings);
    REQUIRE(mcs == Approx(420));
  }
  SECTION("2D")
  {
    double read[4] {2.25, 24.5, 1.225, 101235/0.4 + 0.5*1.225*500};
    double mcs = cartdg::cpg_euler_max<4, 1, 1>(read, 1, settings);
    REQUIRE(mcs == Approx(360).epsilon(0.01));
  }
}

TEST_CASE("minimum ratio of thermodynamic variables")
{
  cartdg::Kernel_settings settings;
  settings.cpg_heat_rat = 1.4;
  settings.max_difference = 0.5;

  SECTION("1D pressure difference")
  {
    double pressure_r [] {1.e5, 1.e5, 0.5e5};
    double pressure_w [] {1.e5, 0.25e5, 0.4e5};
    double state_r [3][3];
    double state_w [3][3];
    for (int i = 0; i < 3; ++i)
    {
      state_r[0][i] = 100.;
      state_w[0][i] = 100.;
      state_r[1][i] = 1.;
      state_w[1][i] = 1.;
      state_r[2][i] = pressure_r[i]/0.4 + 0.5*1.*100*100;
      state_w[2][i] = pressure_w[i]/0.4 + 0.5*1.*100*100;
    }
    double step = cartdg::cpg_euler_physical_step<3, 3, 1>(&state_r[0][0], &state_w[0][0], 1, settings);
    REQUIRE(step == Approx(2./3.));
  }

  SECTION("1D density difference")
  {
    double mass_r [50] {1., 1., 1.};
    double mass_w [50] {0.4, 4., 0.5};
    for (int i = 3; i < 50; ++i) {mass_r[i] = 1.; mass_w[i] = 1.;};
    double state_r [10][3][5];
    double state_w [10][3][5];
    for (int i_elem = 0; i_elem < 10; ++i_elem)
    {
      for (int i = 0; i < 5; ++i)
      {
        state_r[i_elem][0][i] = 100.;
        state_w[i_elem][0][i] = 100.;
        state_r[i_elem][1][i] = mass_r[5*i_elem + i];
        state_w[i_elem][1][i] = mass_w[5*i_elem + i];
        state_r[i_elem][2][i] = 1.e5/0.4 + 0.5*mass_r[5*i_elem + i]*100*100;
        state_w[i_elem][2][i] = 1.e5/0.4 + 0.5*mass_w[5*i_elem + i]*100*100;
      }
    }
    double step = cartdg::cpg_euler_physical_step<3, 5, 1>(&state_r[0][0][0], &state_w[0][0][0], 10, settings);
    REQUIRE(step == Approx(5./6.));
  }

  SECTION("3D large density and small pressure difference")
  {
    double mass_r = 0.5;
    double pressure_r = 2.e5;
    double ener_r = pressure_r/0.4 + 0.5*(100*100 + 200*200 + 50*50)/mass_r;
    double state_r [5] {100, 200, 50, mass_r, ener_r};
    double mass_w = -0.5;
    double pressure_w = 1.5e5;
    double ener_w = pressure_w/0.4 + 0.5*(100*100 + 200*200 + 50*50)/mass_w;
    double state_w [5] {100, 200, 50, mass_w, ener_w};
    double step = cartdg::cpg_euler_physical_step<5, 1, 1>(state_r, state_w, 1, settings);
    REQUIRE(step == Approx(0.25));
  }

  SECTION("3D zero density")
  {
    double mass_r = 0.5;
    double pressure_r = 2.e5;
    double ener_r = pressure_r/0.4 + 0.5*(100*100 + 200*200 + 50*50)/mass_r;
    double state_r [5] {100, 200, 50, mass_r, ener_r};
    double mass_w = 0.;
    double ener_w = 0.9*ener_r;
    double state_w [5] {100, 200, 50, mass_w, ener_w};
    double step = cartdg::cpg_euler_physical_step<5, 1, 1>(state_r, state_w, 1, settings);
    REQUIRE(step <= Approx(0.5));
  }

  SECTION("2D velocity difference")
  {
    double ener = 3*0.5*(100*100 + 200*200)/1.;
    double state_r [4] {100, 200, 1., ener};
    double state_w [4] {200, 400, 1., ener};
    double step = cartdg::cpg_euler_physical_step<4, 1, 1>(state_r, state_w, 1, settings);
    REQUIRE(step <= Approx(0.5));
  }
}