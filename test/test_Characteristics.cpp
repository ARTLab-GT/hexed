#include <catch2/catch_all.hpp>
#include <hexed/pde.hpp>

TEST_CASE("Characteristics")
{
  // arbitrary state to linearize about
  double mass = 1.225;
  hexed::Mat<3> veloc {10., 4., 12.};
  double pres = 101325;
  hexed::Mat<> state(5);
  state(Eigen::seqN(0, 3)) = mass*veloc;
  state(3) = mass;
  state(4) = pres/.4 + .5*mass*veloc.squaredNorm();
  hexed::Mat<3> normal {1., 1., 1.};
  // find eigenvalues
  hexed::pde::Navier_stokes<>::Pde<3>::Characteristics c(state, normal);
  auto eigvals = c.eigvals();
  // decompose another arbitrary state
  hexed::Mat<> state1(5);
  state1 << 1., 2., 3., 1.3, 2e5;
  auto decomp = c.decomp(state1);
  // check that decomposition sums to the input state
  REQUIRE((decomp.rowwise().sum() - state1).cwiseQuotient(state1).norm() == Catch::Approx(0).scale(1.));
  // check that the columns are indeed eigenvectors of linearized flux
  double diff = 1e-6; // use a small perturbation so that flux is effectively linear
  hexed::pde::Navier_stokes<>::Pde<3> ns;
  hexed::pde::Navier_stokes<>::Pde<3>::Computation<1> comp (ns);
  comp.normal = normal;
  for (int i = 0; i < 3; ++i) {
    hexed::Mat<> eigvec = decomp(Eigen::all, i);
    // compute flux perturbation
    hexed::Mat<5> perturb = hexed::Mat<5>::Zero();
    comp.state.setZero();
    comp.state(Eigen::seqN(0, 5)) = state;
    comp.compute_flux_conv();
    perturb -= comp.flux_conv;
    comp.state(Eigen::seqN(0, 5)) += diff*eigvec;
    comp.compute_flux_conv();
    perturb += comp.flux_conv;
    perturb /= normal.norm();
    REQUIRE((perturb.cwiseQuotient(diff*eigvals(i)*eigvec) - hexed::Mat<>::Ones(5)).norm() == Catch::Approx(0.).scale(1.));
  }
}
