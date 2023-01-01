#include <catch2/catch.hpp>
#include <hexed/pde.hpp>

TEST_CASE("Euler")
{
  double mass = 1.225;
  double pressure = 101325;

  SECTION("flux")
  {
    hexed::Mat<2> normal(.6, .8);
    hexed::Mat<2> orth(-.8, .6);
    double nrml_veloc = 7.3;
    double orth_veloc = -4.8;
    hexed::Mat<2> veloc = nrml_veloc*normal + orth_veloc*orth;
    hexed::Mat<4> state {mass*veloc(0), mass*veloc(1), mass, pressure/.4 + .5*veloc.squaredNorm()*mass};
    REQUIRE(hexed::pde::Euler<2>::pressure(state) == Approx(pressure));
    hexed::Mat<4> flux = hexed::pde::Euler<2>::flux(state, normal);
    REQUIRE(flux(0) == Approx(state(0)*nrml_veloc + pressure*normal(0)));
    REQUIRE(flux(1) == Approx(state(1)*nrml_veloc + pressure*normal(1)));
    REQUIRE(flux(2) == Approx(mass*nrml_veloc));
    REQUIRE(flux(3) == Approx(nrml_veloc*(state(3) + pressure)));
  }

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

TEST_CASE("Advection")
{
  SECTION("flux")
  {
    hexed::Mat<2> normal {.6, .8};
    hexed::Mat<2> orth {-.8, .6};
    double nrml_veloc = 7.3;
    double orth_veloc = -4.8;
    hexed::Mat<2> veloc = nrml_veloc*normal + orth_veloc*orth;
    double scalar = 0.81;
    hexed::Mat<3> state {veloc(0), veloc(1), scalar};
    hexed::Mat<1> flux = hexed::pde::Advection<2>::flux(state, normal);
    REQUIRE(flux(0) == Approx(nrml_veloc*scalar));
  }

  SECTION("flux_num")
  {
    auto normal = hexed::Mat<1>::Constant(1.);
    hexed::Mat<2, 2> state;
    state << -3., -2.,
             .21, .23;
    hexed::Mat<1> flux = hexed::pde::Advection<1>::flux_num(state, normal);
    REQUIRE(flux(0) == Approx(-2.*.23*1.));
  }
}
