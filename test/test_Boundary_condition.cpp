#include <catch2/catch.hpp>
#include <cartdgConfig.hpp>
#include <Boundary_condition.hpp>

class Dummy : public cartdg::Boundary_condition
{
  public:
  virtual void apply(cartdg::Boundary_face&) {}
};

TEST_CASE("Typed_boundary_connection")
{
  cartdg::Storage_params params {3, 4, 2, 4};
  cartdg::Element element {params};
  Dummy bc;
  cartdg::Typed_bound_connection<cartdg::Element> tbc0 {element, {1}, bc, false}; // FIXME: address inside_positive for deformed
  cartdg::Typed_bound_connection<cartdg::Element> tbc1 {element, {1}, bc, true};
  REQUIRE(tbc0.n_var() == 4);
  REQUIRE(tbc0.size() == 4*4);
  // check that the correct face of the element is retrieved
  REQUIRE(tbc0.inside_face() == element.face() + (2*1 + 0)*4*4);
  REQUIRE(tbc1.inside_face() == element.face() + (2*1 + 1)*4*4);
  // check that ghost data exists (otherwise segfault)
  tbc0.ghost_face()[0] = 1.;
  tbc0.ghost_face()[4*4 - 1] = 1.;
  // check that order of faces is correct
  REQUIRE(tbc0.face(0) == tbc0.ghost_face());
  REQUIRE(tbc0.face(1) == tbc0.inside_face());
  REQUIRE(tbc1.face(0) == tbc1.inside_face());
  REQUIRE(tbc1.face(1) == tbc1.ghost_face());
}

TEST_CASE("Freestream")
{
  const int row_size = cartdg::config::max_row_size;
  cartdg::Storage_params params {3, 5, 3, row_size};
  cartdg::Element element {params};
  const int n_qpoint = row_size*row_size;
  cartdg::Freestream freestream {{10., 30., -20., 1.3, 1.2e5}};
  cartdg::Typed_bound_connection<cartdg::Element> tbc {element, {1}, freestream, false};
  // set inside face to something arbitrary
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    tbc.inside_face()[0*n_qpoint + i_qpoint] = 20.;
    tbc.inside_face()[1*n_qpoint + i_qpoint] = -10.;
    tbc.inside_face()[2*n_qpoint + i_qpoint] = 50.;
    tbc.inside_face()[3*n_qpoint + i_qpoint] = 0.9;
    tbc.inside_face()[4*n_qpoint + i_qpoint] = 1e4;
  }
  freestream.apply(tbc);
  // check that ghost face is equal to freestream
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    REQUIRE(tbc.ghost_face()[0*n_qpoint + i_qpoint] == Approx(10.));
    REQUIRE(tbc.ghost_face()[1*n_qpoint + i_qpoint] == Approx(30.));
    REQUIRE(tbc.ghost_face()[2*n_qpoint + i_qpoint] == Approx(-20.));
    REQUIRE(tbc.ghost_face()[3*n_qpoint + i_qpoint] == Approx(1.3));
    REQUIRE(tbc.ghost_face()[4*n_qpoint + i_qpoint] == Approx(1.2e5));
  }
}

TEST_CASE("Nonpenetration")
{
  const int row_size = cartdg::config::max_row_size;
  cartdg::Storage_params params {3, 4, 2, row_size};
  cartdg::Deformed_element element {params};
  cartdg::Nonpenetration nonpen;
  cartdg::Typed_bound_connection<cartdg::Deformed_element> tbc {element, {{0, 0}, {1, 0}}, nonpen, true};
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    double qpoint_jacobian [] {1., 3., 0., 4.};
    for (int i_jac = 0; i_jac < 4; ++i_jac) tbc.jacobian()[i_jac*row_size + i_qpoint] = qpoint_jacobian[i_jac];
    double state [] {1., 1., 1.2, 1e5/0.4 + 0.5*1.2*2.};
    for (int i_var = 0; i_var < 4; ++i_var) tbc.inside_face()[(4 + i_var)*row_size + i_qpoint] = state[i_var];
  }
  tbc.inside_face()[4*row_size] = 0.; // momentum[0] at qpoint[0]
  tbc.inside_face()[5*row_size] = 0.; // momentum[1] at qpoint[0]
  nonpen.apply(tbc);
  double normal_momentum_flux = 4./5.*tbc.inside_face()[4*row_size] - 3./5.*tbc.inside_face()[5*row_size];
  double tangential_momentum_flux = 3./5.*tbc.inside_face()[4*row_size] + 4./5.*tbc.inside_face()[5*row_size];
  REQUIRE(normal_momentum_flux == Approx(1e5*5.)); // when momentum == 0, normal momentum flux is pressure
  REQUIRE(tangential_momentum_flux == Approx(0.).scale(1.)); // when momentum != 0, tangential momentum flux is 0
  REQUIRE(tbc.inside_face()[6*row_size + 1] == Approx(0.).scale(1.)); // when momentum != 0, mass flux is 0
}
