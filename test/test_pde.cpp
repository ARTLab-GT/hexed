#include <catch2/catch.hpp>
#include <hexed/pde.hpp>

TEST_CASE("Navier_stokes")
{
  double mass = 1.225;
  double pressure = 101325;

  SECTION("flux")
  {
    hexed::Mat<2> normal {.6, .8};
    hexed::Mat<2> orth {-.8, .6};
    double nrml_veloc = 7.3;
    double orth_veloc = -4.8;
    hexed::Mat<2> veloc = nrml_veloc*normal + orth_veloc*orth;
    normal *= .9;
    hexed::Mat<4> state {mass*veloc(0), mass*veloc(1), mass, pressure/.4 + .5*veloc.squaredNorm()*mass};
    REQUIRE(hexed::pde::Navier_stokes<false>::Pde<2>::pressure(state) == Approx(pressure));
    hexed::Mat<4> flux = hexed::pde::Navier_stokes<false>::Pde<2>::flux(state, normal);
    REQUIRE(flux(0)    == Approx(.9*state(0)*nrml_veloc + pressure*normal(0)));
    REQUIRE(flux(1)    == Approx(.9*state(1)*nrml_veloc + pressure*normal(1)));
    REQUIRE(flux(2)/.9 == Approx(mass*nrml_veloc));
    REQUIRE(flux(3)/.9 == Approx(nrml_veloc*(state(3) + pressure)));
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
      normal *= 1.3;
      flux = hexed::pde::Navier_stokes<false>::Pde<3>::flux_num(state, normal);
      correct = hexed::pde::Navier_stokes<false>::Pde<3>::flux(state(Eigen::all, 0), normal);
      REQUIRE((flux - correct).norm() == Approx(0).scale(1.));

      normal.setUnit(1);
      flux = hexed::pde::Navier_stokes<false>::Pde<3>::flux_num(state, normal);
      correct = hexed::pde::Navier_stokes<false>::Pde<3>::flux(state(Eigen::all, 1), normal);
      REQUIRE((flux - correct).norm() == Approx(0).scale(1.));

      normal.setOnes();
      flux = hexed::pde::Navier_stokes<false>::Pde<3>::flux_num(state, normal);
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
      flux = hexed::pde::Navier_stokes<false>::Pde<3>::flux_num(state, normal);
      REQUIRE(flux(3) == Approx(0.).scale(1.));
      REQUIRE(flux(4) == Approx(0.).scale(1.));
    }
  }

  SECTION("flux_viscous")
  {
    hexed::Mat<4> vrhot; // contains velocity, density, and temperature in order
    hexed::Mat<2, 4> vrhot_grad; // gradients of above
    vrhot << 10., 11., 1.2, 300.;
    vrhot_grad << .2, -.3, .1, 10.,
                  -.1, .5, .01, 13;
    auto state = vrhot; // conserved variables
    auto seq = Eigen::seqN(0, 2);
    state(seq) *= vrhot(2);
    double cv = 287.058/.4;
    state(3) = vrhot(2)*cv*vrhot(3) + .5*state(seq).dot(vrhot(seq));
    // compute gradient of conserved variables
    hexed::Mat<2, 4> state_grad;
    for (int i = 0; i < 2; ++i) {
      state_grad(Eigen::all, i) = vrhot(2)*vrhot_grad(Eigen::all, i) + vrhot(i)*vrhot_grad(Eigen::all, 2);
    }
    state_grad(Eigen::all, 2) = vrhot_grad(Eigen::all, 2);
    state_grad(Eigen::all, 3) = cv*(vrhot_grad(Eigen::all, 2)*vrhot(3) + vrhot_grad(Eigen::all, 3)*vrhot(2))
                                + .5*(state_grad(Eigen::all, seq)*vrhot(seq) + vrhot_grad(Eigen::all, seq)*state(seq));
    // now compute the correct viscous flux (with the original variables to help)
    // and check that `Navier_stokes` gets the same thing
    double visc = 1.81206e-5;
    double cond = visc*1.4*cv/.71;
    hexed::Mat<2, 4> correct = -.9*state_grad;
    auto veloc_grad = vrhot_grad(Eigen::all, seq);
    hexed::Mat<2, 2> stress = visc*(veloc_grad + veloc_grad.transpose() - 2./3.*veloc_grad.trace()*hexed::Mat<2, 2>::Identity());
    correct(Eigen::all, seq) -= stress;
    correct(Eigen::all, 3) -= stress*vrhot(seq) + cond*vrhot_grad(Eigen::all, 3);
    REQUIRE((hexed::pde::Navier_stokes<true>::Pde<2>::flux_visc(state, state_grad, .9) - correct).norm() < 1e-10);
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
    normal *= .5;
    double scalar = 0.81;
    hexed::Mat<3> state {veloc(0), veloc(1), scalar};
    hexed::Mat<1> flux = hexed::pde::Advection<2>::flux(state, normal);
    REQUIRE(flux(0) == Approx(.5*nrml_veloc*scalar));
  }

  SECTION("flux_num")
  {
    auto normal = hexed::Mat<1>::Constant(.6);
    hexed::Mat<2, 2> state;
    state << -3., -2.,
             .21, .23;
    hexed::Mat<1> flux = hexed::pde::Advection<1>::flux_num(state, normal);
    REQUIRE(flux(0) == Approx(-2.*.23*.6));
  }
}
