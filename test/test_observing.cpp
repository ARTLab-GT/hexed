#include <catch.hpp>

#include <cartdgConfig.hpp>
#include <Kernel_settings.hpp>
#include <get_mcs_cpg_euler.hpp>
#include <Gauss_lobatto.hpp>
#include <observing/indicator.hpp>

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

TEST_CASE("Discontinuity indicator")
{
  const int row_size = MAX_BASIS_ROW_SIZE;
  cartdg::Gauss_lobatto basis (row_size);
  cartdg::Gauss_lobatto basis1 (row_size - 1);

  SECTION("1D")
  {
    SECTION("different types of functions")
    {
      double read [row_size];
      auto weights = basis.node_weights();
      auto ortho = basis.orthogonal(row_size - 1);

      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
      {
        read[i_qpoint] = std::exp(basis.node(i_qpoint)*0.5);
      }
      REQUIRE(cartdg::indicator<row_size, row_size>(read, weights.data(), ortho.data()) == 0.);
      read[1] = 2.;
      REQUIRE(cartdg::indicator<row_size, row_size>(read, weights.data(), ortho.data()) == 1.);

      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
      {
        read[i_qpoint] = 1;
      }
      REQUIRE(cartdg::indicator<row_size, row_size>(read, weights.data(), ortho.data()) == 0.);

      for (int i_qpoint = 2; i_qpoint < row_size; ++i_qpoint)
      {
        read[i_qpoint] = 1.5;
      }
      REQUIRE(cartdg::indicator<row_size, row_size>(read, weights.data(), ortho.data()) == 1.);

      SECTION("smoothness of indicator function")
      {
        for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
        {
          read[i_qpoint] = 1.;
        }
        double prev = 0.;
        double curr;
        for (int pow = -200; pow < 0; ++pow)
        {
          read[0] = 1. + std::exp(pow*0.05);
          curr = cartdg::indicator<row_size, row_size>(read, weights.data(), ortho.data());
          REQUIRE(curr - prev < 0.1);
          REQUIRE(curr - prev >= 0.);
          prev = curr;
        }
        REQUIRE(curr == 1.);
      }
    }
    SECTION("even/odd")
    {
      double read [row_size - 1];
      auto weights = basis1.node_weights();
      auto ortho = basis1.orthogonal(row_size - 2);
      for (int i_qpoint = 0; i_qpoint < row_size - 1; ++i_qpoint)
      {
        read[i_qpoint] = std::exp(basis1.node(i_qpoint)*0.5);
      }
      REQUIRE(cartdg::indicator<row_size - 1, row_size - 1>(read, weights.data(), ortho.data()) == 0.);
      read[1] = 2.;
      REQUIRE(cartdg::indicator<row_size - 1, row_size - 1>(read, weights.data(), ortho.data()) == 1.);
    }
  }

  SECTION("3D")
  {
    double read [row_size][row_size][row_size];
    auto weights = basis.node_weights();
    auto ortho = basis.orthogonal(row_size - 1);

    const int n_qpoint = row_size*row_size*row_size;
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      (&read[0][0][0])[i_qpoint] = 1.;
    }
    REQUIRE(cartdg::indicator<n_qpoint, row_size>(read[0][0], weights.data(), ortho.data()) == 0.);

    SECTION("all directions activated")
    {
      read[1][2][1] = 2.;
      REQUIRE(cartdg::indicator<n_qpoint, row_size>(read[0][0], weights.data(), ortho.data()) == 1.);
    }

    SECTION("one direction activated")
    {
      for (int i_row = 0; i_row < row_size; ++i_row)
      {
        for (int j_row = 0; j_row < row_size; ++j_row)
        {
          read[i_row][1][j_row] = 0.9;
        }
      }
      REQUIRE(cartdg::indicator<n_qpoint, row_size>(read[0][0], weights.data(), ortho.data()) == 1.);
    }
  }
}
