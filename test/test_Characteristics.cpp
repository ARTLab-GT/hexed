#include <iostream>
#include <catch2/catch.hpp>
#include <hexed/Characteristics.hpp>
#include <hexed/pde.hpp>

TEST_CASE("Characteristics")
{
  double mass = 1.225;
  //hexed::Mat<3> veloc {10., 4., 12.};
  hexed::Mat<3> veloc {10., 10., 10.};
  double pres = 101325;
  hexed::Mat<> state(5);
  state(Eigen::seqN(0, 3)) = mass*veloc;
  state(3) = mass;
  state(4) = pres/.4 + .5*mass*veloc.squaredNorm();
  hexed::Mat<3> normal {1., 1., 1.};
  hexed::Characteristics c(state, normal);
  auto eigvals = c.eigvals();
  hexed::Mat<> state1(5);
  state1 << 1., 2., 3., 1.3, 2e5;
  auto decomp = c.decomp(state1);
  REQUIRE(((decomp.rowwise().sum() - state1).cwiseQuotient(state1)).norm() == Approx(0).scale(1.));
  double diff = 1e-7;
  hexed::pde::Navier_stokes<>::Pde<3> ns;
  for (int i = 0; i < 3; ++i) {
    printf("%i\n", i);
    hexed::Mat<> eigvec = decomp(Eigen::all, i);
    hexed::Mat<> perturb = (ns.flux(state + diff*eigvec, normal) - ns.flux(state, normal))/normal.norm();
    CHECK((perturb.cwiseQuotient(diff*eigvals(i)*eigvec) - hexed::Mat<>::Ones(5)).norm() == Approx(0.).scale(1.));
  }
}
