#include <catch2/catch_all.hpp>
#include <hexed/config.hpp>
#include <hexed/Boundary_condition.hpp>
#include <hexed/Spacetime_func.hpp>
#include <hexed/constants.hpp>
#include <hexed/connection.hpp>

class Dummy : public hexed::Boundary_condition
{
  public:
  virtual void apply_state(hexed::Boundary_face&) {}
};

TEST_CASE("Typed_boundary_connection")
{
  hexed::Storage_params params {3, 4, 2, 4};
  hexed::Element element {params};
  Dummy bc;
  hexed::Typed_bound_connection<hexed::Element> tbc0 {element, 1, false, 0};
  REQUIRE(element.faces[2] == tbc0.state());
  REQUIRE(tbc0.ghost_face() == tbc0.state() + params.n_dof()/params.row_size);
  hexed::Typed_bound_connection<hexed::Element> tbc1 {element, 1,  true, 1};
  REQUIRE(element.faces[3] == tbc1.state());
  REQUIRE(tbc1.ghost_face() == tbc1.state() + params.n_dof()/params.row_size);
  REQUIRE(tbc0.storage_params().n_var == 4);
  // check that the correct face of the element is retrieved
  REQUIRE(tbc0.inside_face() == element.faces[2*1 + 0]);
  REQUIRE(tbc1.inside_face() == element.faces[2*1 + 1]);
  // check that ghost data exists (otherwise segfault)
  tbc0.ghost_face()[0] = 1.;
  tbc0.ghost_face()[4*4 - 1] = 1.;
  // check that order of faces is correct
  REQUIRE(tbc0.state()                                  == tbc0.inside_face());
  REQUIRE(tbc0.state() + params.n_dof()/params.row_size == tbc0.ghost_face());
  REQUIRE(tbc1.state()                                  == tbc1.inside_face());
  REQUIRE(tbc1.state() + params.n_dof()/params.row_size == tbc1.ghost_face());
  // check that the direction info is correct
  REQUIRE(tbc0.direction().i_dim[0] == 1);
  REQUIRE(tbc0.direction().i_dim[1] == 1);
  REQUIRE(tbc0.direction().face_sign[0] == 0);
  REQUIRE(tbc0.direction().face_sign[1] == 1);
  REQUIRE(tbc1.direction().face_sign[0] == 1);
  REQUIRE(tbc1.direction().face_sign[1] == 0);
}

TEST_CASE("Freestream")
{
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 5, 3, row_size};
  hexed::Element element {params};
  const int n_qpoint = row_size*row_size;
  hexed::Freestream freestream {hexed::Mat<5>{10., 30., -20., 1.3, 1.2e5}};
  hexed::Typed_bound_connection<hexed::Element> tbc {element, 1, false, 0};
  // set inside face to something arbitrary
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    tbc.inside_face()[0*n_qpoint + i_qpoint] = 20.;
    tbc.inside_face()[1*n_qpoint + i_qpoint] = -10.;
    tbc.inside_face()[2*n_qpoint + i_qpoint] = 50.;
    tbc.inside_face()[3*n_qpoint + i_qpoint] = 0.9;
    tbc.inside_face()[4*n_qpoint + i_qpoint] = 1e4;
  }
  freestream.apply_state(tbc);
  // check that ghost face is equal to freestream
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    REQUIRE(tbc.ghost_face()[0*n_qpoint + i_qpoint] == Catch::Approx(10.));
    REQUIRE(tbc.ghost_face()[1*n_qpoint + i_qpoint] == Catch::Approx(30.));
    REQUIRE(tbc.ghost_face()[2*n_qpoint + i_qpoint] == Catch::Approx(-20.));
    REQUIRE(tbc.ghost_face()[3*n_qpoint + i_qpoint] == Catch::Approx(1.3));
    REQUIRE(tbc.ghost_face()[4*n_qpoint + i_qpoint] == Catch::Approx(1.2e5));
  }
}

TEST_CASE("Riemann_invariants")
{
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 5, 3, row_size};
  hexed::Element element {params};
  const int n_qpoint = row_size*row_size;
  hexed::Mat<5> fs {10., 30., -20., 1.3, 4e5};
  hexed::Riemann_invariants ri {fs};
  hexed::Typed_bound_connection<hexed::Element> tbc {element, 1, false, 0};
  hexed::Mat<5> inside_state {1/1.2, -600/1.2, 1/1.2, 1.2, 101325/.4 + .5*1.2*360002};
  SECTION("supersonic inflow")
  {
    // set the first point to supersonic outflow and the rest to supersonic inflow
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      hexed::Mat<5> state = inside_state;
      if (i_qpoint) state(1) *= -1;
      for (int i_var = 0; i_var < 5; ++i_var) {
        tbc.inside_face()[i_var*n_qpoint + i_qpoint] = state(i_var);
      }
      for (int i_dim = 0; i_dim < 3; ++i_dim) tbc.surface_normal()[i_dim*n_qpoint + i_qpoint] = (i_dim == 1);
    }
    ri.apply_state(tbc);
    // test that the first point is the inside state and the rest are the freestream state
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        REQUIRE(tbc.ghost_face()[i_var*n_qpoint + i_qpoint] == Catch::Approx(i_qpoint ? fs[i_var] : inside_state(i_var)));
      }
    }
    // set an arbitrary viscous flux
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        tbc.inside_face()[(10 + i_var)*n_qpoint + i_qpoint] = 1.;
      }
    }
    ri.apply_flux(tbc);
    // test that the flux is left alone at the points where the state was set and set to zero where it was not
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        REQUIRE(tbc.ghost_face()[(10 + i_var)*n_qpoint + i_qpoint] == Catch::Approx(i_qpoint ? 1. : 0.).scale(1.));
      }
    }
  }
}

TEST_CASE("Function_bc")
{
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {2, 4, 2, row_size};
  hexed::Element element {params};
  const int n_qpoint = row_size;
  hexed::Annular_diffusion_test func(1.7, 2., 1e5);
  hexed::Function_bc bc(func);
  hexed::Typed_bound_connection<hexed::Element> tbc {element, 1, false, 0};
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    // set inside face to something arbitrary
    tbc.inside_face()[0*n_qpoint + i_qpoint] = 20.;
    tbc.inside_face()[1*n_qpoint + i_qpoint] = -10.;
    tbc.inside_face()[2*n_qpoint + i_qpoint] = 5.;
    tbc.inside_face()[3*n_qpoint + i_qpoint] = 1e4;
  }
  // set position at node 0 to have radius 2*e (not the actual position, but for this test that doesn't matter)
  tbc.surface_position()[0*n_qpoint + 0] =  1.2*std::exp(1.);
  tbc.surface_position()[1*n_qpoint + 0] = -1.6*std::exp(1.);
  // set position at node 4 to have radius 2*e^2
  tbc.surface_position()[0*n_qpoint + 4] =  1.2*std::exp(2.);
  tbc.surface_position()[1*n_qpoint + 4] =  1.6*std::exp(2.);
  bc.apply_state(tbc);
  // check that ghost face state is correct at the qpoints where position was set
  REQUIRE(tbc.ghost_face()[0*n_qpoint + 0] == Catch::Approx(0.).scale(1.));
  REQUIRE(tbc.ghost_face()[1*n_qpoint + 0] == Catch::Approx(0.).scale(1.));
  REQUIRE(tbc.ghost_face()[2*n_qpoint + 0] == Catch::Approx(1.7));
  REQUIRE(tbc.ghost_face()[3*n_qpoint + 0] == Catch::Approx(1e5));
  REQUIRE(tbc.ghost_face()[2*n_qpoint + 4] == Catch::Approx(3.4));
}

TEST_CASE("Nonpenetration")
{
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 4, 2, row_size};
  hexed::Deformed_element element {params};
  hexed::Nonpenetration nonpen;
  hexed::Typed_bound_connection<hexed::Deformed_element> tbc {element, 0, true, 0};
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    double qpoint_nrml [] {-4., 3.};
    for (int i_dim = 0; i_dim < 2; ++i_dim) tbc.normal()[i_dim*row_size + i_qpoint] = qpoint_nrml[i_dim];
    double state [] {1., 1., 1.2, 1e5/0.4 + 0.5*1.2*2.};
    for (int i = 0; i < 2; ++i) {
      for (int i_var = 0; i_var < 4; ++i_var) tbc.inside_face()[(8*i + i_var)*row_size + i_qpoint] = state[i_var];
    }
  }
  SECTION("apply_state")
  {
    nonpen.apply_state(tbc);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      // require tangential momentum unchanged
      REQUIRE(  3*tbc.ghost_face()[0*row_size + i_qpoint]
              + 4*tbc.ghost_face()[1*row_size + i_qpoint] == Catch::Approx(7.));
      // require normal momentum flipped
      REQUIRE( -4*tbc.ghost_face()[0*row_size + i_qpoint]
              + 3*tbc.ghost_face()[1*row_size + i_qpoint] == Catch::Approx(1.));
    }
  }
  SECTION("apply_flux")
  {
    nonpen.apply_flux(tbc);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      // require tangential momentum flux flipped
      REQUIRE(  3*tbc.ghost_face()[8*row_size + i_qpoint]
              + 4*tbc.ghost_face()[9*row_size + i_qpoint] == Catch::Approx(-7.));
      // require normal momentum flux unchanged
      REQUIRE( -4*tbc.ghost_face()[8*row_size + i_qpoint]
              + 3*tbc.ghost_face()[9*row_size + i_qpoint] == Catch::Approx(-1.));
      // require scalar flux flipped
      REQUIRE(tbc.ghost_face()[10*row_size + i_qpoint] == Catch::Approx(-tbc.inside_face()[10*row_size + i_qpoint]));
      REQUIRE(tbc.ghost_face()[11*row_size + i_qpoint] == Catch::Approx(-tbc.inside_face()[11*row_size + i_qpoint]));
    }
  }
}

TEST_CASE("No_slip")
{
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 4, 2, row_size};
  hexed::Deformed_element element {params};
  hexed::Typed_bound_connection<hexed::Deformed_element> tbc {element, 0, false, 0};
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    for (int i_dim = 0; i_dim < 2; ++i_dim) tbc.normal()[i_dim*row_size + i_qpoint] = .7/std::sqrt(2.);
  }
  double state [] {1., 1., 1.2, 1e5/0.4 + 0.5*1.2*2.};
  double flux [] {10., -20., 1.3, 10.};
  SECTION("isothermal")
  {
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) tbc.inside_face()[i_var*row_size + i_qpoint] = state[i_var];
    }
    hexed::No_slip no_slip(hexed::No_slip::internal_energy, 1e6);
    no_slip.apply_state(tbc);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_dim = 0; i_dim < 2; ++i_dim) REQUIRE(tbc.ghost_face()[i_dim*row_size + i_qpoint] == Catch::Approx(-1.));
      REQUIRE(tbc.ghost_face()[2*row_size + i_qpoint] == Catch::Approx(1.2));
      REQUIRE(std::sqrt(tbc.ghost_face()[3*row_size + i_qpoint]*tbc.inside_face()[3*row_size + i_qpoint]) == Catch::Approx(1e6*1.2));
    }
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) tbc.inside_face()[(8 + i_var)*row_size + i_qpoint] = flux[i_var];
    }
    no_slip.apply_flux(tbc);
    flux[2] *= -1; // mass flux should always inverted
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(tbc.ghost_face()[(8 + i_var)*row_size + i_qpoint] == Catch::Approx(flux[i_var]));
      }
    }
  }
  SECTION("specified flux")
  {
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) tbc.inside_face()[i_var*row_size + i_qpoint] = state[i_var];
    }
    hexed::No_slip no_slip(hexed::No_slip::heat_flux, 3.);
    no_slip.apply_state(tbc);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_dim = 0; i_dim < 2; ++i_dim) REQUIRE(tbc.ghost_face()[i_dim*row_size + i_qpoint] == Catch::Approx(-1.));
      REQUIRE(tbc.ghost_face()[2*row_size + i_qpoint] == Catch::Approx(state[2]));
      REQUIRE(tbc.ghost_face()[3*row_size + i_qpoint] == Catch::Approx(state[3]));
    }
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) tbc.inside_face()[(8 + i_var)*row_size + i_qpoint] = flux[i_var];
    }
    no_slip.apply_flux(tbc);
    flux[2] *= -1; // mass flux should always inverted
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 3; ++i_var) {
        REQUIRE(tbc.ghost_face()[(8 + i_var)*row_size + i_qpoint] == Catch::Approx(flux[i_var]));
      }
      REQUIRE((tbc.ghost_face()[11*row_size + i_qpoint] + tbc.inside_face()[11*row_size + i_qpoint])/2 == Catch::Approx(-3.*.7));
    }
  }
  SECTION("specified emissivity")
  {
    state[3] = 1e5/.4;
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) {
        tbc.inside_face()[i_var*row_size + i_qpoint] = state[i_var];
        tbc.inside_face()[(8 + i_var)*row_size + i_qpoint] = flux[i_var];
      }
    }
    hexed::No_slip no_slip(hexed::No_slip::emissivity, .8);
    double temp = 1e5/1.2/hexed::constants::specific_gas_air;
    no_slip.apply_state(tbc);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) {
        tbc.inside_face()[i_var*row_size + i_qpoint] = 0;
      }
    }
    no_slip.apply_flux(tbc);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      REQUIRE((tbc.ghost_face()[11*row_size + i_qpoint] + tbc.inside_face()[11*row_size + i_qpoint])/2 == Catch::Approx(-.8*hexed::constants::stefan_boltzmann*std::pow(temp, 4)*.7));
    }
  }
}

TEST_CASE("Copy")
{
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 5, 3, row_size};
  hexed::Element element {params};
  const int n_qpoint = row_size*row_size;
  hexed::Copy copy;
  hexed::Typed_bound_connection<hexed::Element> tbc {element, 1, false, 0};
  // set inside face to something arbitrary
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    tbc.inside_face()[0*n_qpoint + i_qpoint] = 20.;
    tbc.inside_face()[1*n_qpoint + i_qpoint] = -10.;
    tbc.inside_face()[2*n_qpoint + i_qpoint] = 50.;
    tbc.inside_face()[3*n_qpoint + i_qpoint] = 0.9;
    tbc.inside_face()[4*n_qpoint + i_qpoint] = 1e4;
  }
  copy.apply_state(tbc);
  // check that ghost face is equal to inside
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    REQUIRE(tbc.ghost_face()[0*n_qpoint + i_qpoint] == Catch::Approx(20.));
    REQUIRE(tbc.ghost_face()[1*n_qpoint + i_qpoint] == Catch::Approx(-10.));
    REQUIRE(tbc.ghost_face()[2*n_qpoint + i_qpoint] == Catch::Approx(50.));
    REQUIRE(tbc.ghost_face()[3*n_qpoint + i_qpoint] == Catch::Approx(0.9));
    REQUIRE(tbc.ghost_face()[4*n_qpoint + i_qpoint] == Catch::Approx(1e4));
  }
}
