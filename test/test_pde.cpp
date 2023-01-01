#include <catch2/catch.hpp>
#include <hexed/pde.hpp>

TEST_CASE("Euler")
{
  double mass = 1.225;
  double pressure = 101325;

  SECTION("flux_num")
  {
    hexed::Mat<5, 2> state;
    hexed::Mat<3> normal;
    hexed::Mat<5> flux;
    hexed::Mat<5> correct;
    SECTION("Reasonable flow")
    {
      hexed::Mat<5, 2> state;
      hexed::Mat<3> normal;
      hexed::Mat<5> flux;
      hexed::Mat<5> correct;
      double velocity0 [] {3*340, 2*340};
      double velocity1 [] {-2*340, -3*340};
      for (int i_side = 0; i_side < 2; ++i_side) {
        state(0, i_side) = mass*velocity0[i_side];
        state(1, i_side) = mass*velocity1[i_side];
        state(2, i_side) = 0;
        state(3, i_side) = mass;
        state(4, i_side) = pressure/0.4 + 0.5*mass*(  velocity0[i_side]*velocity0[i_side]
                                                    + velocity1[i_side]*velocity1[i_side]);
      }
      normal.setUnit(0);
      flux = hexed::pde::Euler<3>::flux_num(state, normal);
      correct = hexed::pde::Euler<3>::flux(state(Eigen::all, 0), normal);
      REQUIRE((flux - correct).norm() == Approx(0).scale(1.));

      normal.setUnit(1);
      flux = hexed::pde::Euler<3>::flux_num(state, normal);
      correct = hexed::pde::Euler<3>::flux(state(Eigen::all, 1), normal);
      REQUIRE((flux - correct).norm() == Approx(0).scale(1.));

      normal.setOnes(1);
      flux = hexed::pde::Euler<3>::flux_num(state, normal);
      REQUIRE(flux(3) == Approx(0.).scale(1.));
    }

    SECTION("Opposing supersonic flows")
    {
      double velocity0 [] {680, -680};
      for (int i_side = 0; i_side < 2; ++i_side) {
        state(0, i_side) = mass*velocity0[i_side];
        state(1, i_side) = 0;
        state(2, i_side) = 0;
        state(3, i_side) = mass;
        state(4, i_side) = pressure/0.4 + 0.5*mass*velocity0[i_side]*velocity0[i_side];
      }
      normal.setUnit(0);
      flux = hexed::pde::Euler<3>::flux_num(state, normal);
      REQUIRE(flux(3) == Approx(0.).scale(1.));
      REQUIRE(flux(4) == Approx(0.).scale(1.));
    }
  }
}

#if 0
TEST_CASE("hll::advection")
{
  double state [2][4][3];
  double nrml [2][3];
  double scalar [2] {1.3, 0.7};
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (int i_qpoint = 0; i_qpoint < 3; ++i_qpoint) {
      // set surface normal
      nrml[0][i_qpoint] = .4;
      nrml[1][i_qpoint] = .5;
      // set advection velocity
      state[i_side][0][i_qpoint] = -.1;
      state[i_side][1][i_qpoint] = -.9;
      // set scalar state;
      state[i_side][2][i_qpoint] = scalar[i_side];
    }
  }
  hexed::hll::advection<2, 3>(state[0][0], nrml[0]);
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (int i_qpoint = 0; i_qpoint < 3; ++i_qpoint) {
      REQUIRE(state[i_side][0][i_qpoint] == Approx(-.1));
      REQUIRE(state[i_side][1][i_qpoint] == Approx(-.9));
      REQUIRE(state[i_side][2][i_qpoint] == Approx((-.1*.4 - .9*.5)*0.7));
    }
  }
}
#endif
