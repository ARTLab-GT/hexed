#include <catch2/catch_all.hpp>
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
    REQUIRE(hexed::pde::Navier_stokes<false>::Pde<2>().pressure(state) == Catch::Approx(pressure));
    hexed::Mat<4> flux = hexed::pde::Navier_stokes<false>::Pde<2>().flux(state, normal);
    REQUIRE(flux(0)    == Catch::Approx(.9*state(0)*nrml_veloc + pressure*normal(0)));
    REQUIRE(flux(1)    == Catch::Approx(.9*state(1)*nrml_veloc + pressure*normal(1)));
    REQUIRE(flux(2)/.9 == Catch::Approx(mass*nrml_veloc));
    REQUIRE(flux(3)/.9 == Catch::Approx(nrml_veloc*(state(3) + pressure)));
  }

  SECTION("flux_viscous")
  {
    // this isn't a particularly good test,
    // since the test and the actual function share a fair amount of code in common.
    // However, i can't think of a really effective way to test this,
    // and at least the gradients start out in terms of different variables
    hexed::Mat<4> vrhot; // contains velocity, density, and temperature in order
    hexed::Mat<2, 4> vrhot_grad; // gradients of above
    vrhot << 10., 11., 1.2, 300.;
    vrhot_grad << .2, -.3, .1, 10.,
                  -.1, .5, .01, 13;
    auto state = vrhot; // conserved variables
    auto seq = Eigen::seqN(0, 2);
    state(seq) *= vrhot(2);
    double cv = hexed::constants::specific_gas_air/.4;
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
    hexed::Mat<2, 4> correct;
    correct.setZero();
    auto veloc_grad = vrhot_grad(Eigen::all, seq);
    hexed::Mat<2, 2> stress = .9*vrhot(2)*veloc_grad.trace()*hexed::Mat<2, 2>::Identity();
    correct(Eigen::all, seq) = -stress;
    correct(Eigen::all, 3) = -stress*vrhot(seq);
    REQUIRE((hexed::pde::Navier_stokes<true>::Pde<2>().flux_visc(state, state_grad, .9) - correct).norm() < 1e-10);
    stress += visc*(veloc_grad + veloc_grad.transpose() - 2./3.*veloc_grad.trace()*hexed::Mat<2, 2>::Identity());
    correct(Eigen::all, seq) = -stress;
    correct(Eigen::all, 3) = -stress*vrhot(seq) - cond*vrhot_grad(Eigen::all, 3);
    REQUIRE((hexed::pde::Navier_stokes<true>::Pde<2>(hexed::Transport_model::constant(visc), hexed::Transport_model::constant(cond)).flux_visc(state, state_grad, .9) - correct).norm() < 1e-10);
  }
}
