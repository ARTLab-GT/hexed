#include <catch.hpp>

#include <kernels/neighbor/cpg_euler_fbc.hpp>

TEST_CASE("flux bc")
{
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 0.1;

  double veloc0 = 2.;
  double veloc1 = 1.;
  double veloc2 = 0.;
  double mass = 1.2;
  double pressure = 1e5;
  double ener = (pressure/(settings.cpg_heat_rat - 1.) + mass*0.5*14);

  double read [2][5][8];
  double write [2][5][8] {};
  double jacobian [2][3][3][8] {};
  int i_elem [] {1, 0, 1};
  int i_dim [] {1, 2, 1};
  int is_positive [] {0, 1, 1};

  for (int i_elem : {0, 1})
  {
    for (int i_qpoint = 0; i_qpoint < 8; ++i_qpoint)
    {
      read[i_elem][0][i_qpoint] = veloc0*mass;
      read[i_elem][1][i_qpoint] = veloc1*mass;
      read[i_elem][3][i_qpoint] = veloc2*mass;
      read[i_elem][3][i_qpoint] = mass;
      read[i_elem][4][i_qpoint] = ener;

      jacobian[i_elem][2][0][i_qpoint] = 1;
      jacobian[i_elem][0][1][i_qpoint] = 1;
      jacobian[i_elem][1][2][i_qpoint] = 1;
      jacobian[i_elem][2][1][i_qpoint] = 2;
    }
  }

  cartdg::cpg_euler_fbc<5, 8, 2>(&read[0][0][0], &write[0][0][0], &jacobian[0][0][0][0],
                                 i_elem, i_dim, is_positive, 2, 2., settings);

  for (int i_qpoint = 0; i_qpoint < 8; ++i_qpoint)
  {
    REQUIRE(write[0][0][i_qpoint] == Approx(0.));
    REQUIRE(write[0][1][i_qpoint] == Approx(0.));
    REQUIRE(write[0][3][i_qpoint] == Approx(0.));
    REQUIRE(write[0][3][i_qpoint] == Approx(0.));
    REQUIRE(write[0][4][i_qpoint] == Approx(0.));
  }

  double veloc_normal = -3;
  double mult = settings.d_t_by_d_pos/2.;
  for (int i_qpoint : {0, 1, 4, 5})
  {
    REQUIRE(write[1][3][i_qpoint] == Approx(-veloc_normal*mass*mult));
    REQUIRE(write[1][4][i_qpoint] == Approx(-veloc_normal*(ener + pressure)*mult));
  }
  for (int i_qpoint : {2, 3, 6, 7})
  {
    REQUIRE(write[1][3][i_qpoint] == Approx(veloc_normal*mass*mult));
    REQUIRE(write[1][4][i_qpoint] == Approx(veloc_normal*(ener + pressure)*mult));
  }
}