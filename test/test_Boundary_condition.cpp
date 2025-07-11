#include <catch2/catch_all.hpp>
#include <hexed/config.hpp>
#include <hexed/Boundary_condition.hpp>
#include <hexed/Spacetime_func.hpp>
#include <hexed/constants.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Simplex_geom.hpp>

TEST_CASE("Freestream") {
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 5, 3, row_size};
  hexed::Tree tree(3, 1.);
  hexed::Element element {params, tree};
  const int n_qpoint = row_size*row_size;
  hexed::Freestream freestream {hexed::Mat<5>{10., 30., -20., 1.3, 1.2e5}};
  hexed::Boundary_connection con(element.face(2), 0, 0);
  // set inside face to something arbitrary
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    con.inside().flow_state()(0)(0)[i_qpoint] = 20.;
    con.inside().flow_state()(0)(1)[i_qpoint] = -10.;
    con.inside().flow_state()(0)(2)[i_qpoint] = 50.;
    con.inside().flow_state()(0)(3)[i_qpoint] = 0.9;
    con.inside().flow_state()(0)(4)[i_qpoint] = 1e4;
  }
  freestream.apply_state(con);
  // check that ghost face is equal to freestream
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    REQUIRE(con.ghost().flow_state()(0)(0)[i_qpoint] == Catch::Approx(10.));
    REQUIRE(con.ghost().flow_state()(0)(1)[i_qpoint] == Catch::Approx(30.));
    REQUIRE(con.ghost().flow_state()(0)(2)[i_qpoint] == Catch::Approx(-20.));
    REQUIRE(con.ghost().flow_state()(0)(3)[i_qpoint] == Catch::Approx(1.3));
    REQUIRE(con.ghost().flow_state()(0)(4)[i_qpoint] == Catch::Approx(1.2e5));
  }
}

TEST_CASE("Riemann_invariants") {
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 5, 3, row_size};
  hexed::Tree tree(3, 1.);
  hexed::Element element {params, tree};
  const int n_qpoint = row_size*row_size;
  hexed::Mat<5> fs {10., 30., -20., 1.3, 4e5};
  hexed::Riemann_invariants ri {fs};
  hexed::Boundary_connection con(element.face(2), 0, 0);
  hexed::Mat<5> inside_state {1/1.2, -600/1.2, 1/1.2, 1.2, 101325/.4 + .5*1.2*360002};
  SECTION("supersonic inflow")
  {
    // set the first point to supersonic outflow and the rest to supersonic inflow
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      hexed::Mat<5> state = inside_state;
      if (i_qpoint) state(1) *= -1;
      for (int i_var = 0; i_var < 5; ++i_var) {
        con.inside().flow_state()(0)(i_var)[i_qpoint] = state(i_var);
      }
      for (int i_dim = 0; i_dim < 3; ++i_dim) con.normal()(i_dim)[i_qpoint] = (i_dim == 1);
    }
    ri.apply_state(con);
    // test that the first point is the inside state and the rest are the freestream state
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        REQUIRE(con.ghost().flow_state()()(0)(i_var)[i_qpoint] == Catch::Approx(i_qpoint ? fs[i_var] : inside_state(i_var)));
      }
    }
    // set an arbitrary viscous flux
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        con.inside().flow_state()(1)(i_var)[i_qpoint] = 1.;
      }
    }
    ri.apply_flux(con);
    // test that the flux is left alone at the points where the state was set and set to zero where it was not
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        REQUIRE(con.ghost().flow_state()(1)(i_var)[i_qpoint] == Catch::Approx(i_qpoint ? 1. : 0.).scale(1.));
      }
    }
  }
}

TEST_CASE("Function_bc") {
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {2, 4, 2, row_size};
  hexed::Tree tree(2, 1.);
  hexed::Element element {params, tree};
  const int n_qpoint = row_size;
  hexed::Annular_diffusion_test func(1.7, 2., 1e5);
  hexed::Function_bc bc(func);
  hexed::Boundary_connection con(element.face(2), 0, 0);
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    // set inside face to something arbitrary
    con.inside().flow_state()(0)(0)[i_qpoint] = 20.;
    con.inside().flow_state()(0)(1)[i_qpoint] = -10.;
    con.inside().flow_state()(0)(2)[i_qpoint] = 5.;
    con.inside().flow_state()(0)(3)[i_qpoint] = 1e4;
  }
  con.normal() = 0.;
  // set position at node 0 to have radius 2*e (not the actual position, but for this test that doesn't matter)
  con.position()(0)[0] =  1.2*std::exp(1.);
  con.position()(1)[0] = -1.6*std::exp(1.);
  // set position at node 4 to have radius 2*e^2
  con.position()(0)[4] =  1.2*std::exp(2.);
  con.position()(1)[4] =  1.6*std::exp(2.);
  bc.apply_state(con);
  // check that ghost face state is correct at the qpoints where position was set
  REQUIRE(con.ghost().flow_state()(0)(0)[0] == Catch::Approx(0.).scale(1.));
  REQUIRE(con.ghost().flow_state()(0)(1)[0] == Catch::Approx(0.).scale(1.));
  REQUIRE(con.ghost().flow_state()(0)(2)[0] == Catch::Approx(1.7));
  REQUIRE(con.ghost().flow_state()(0)(3)[0] == Catch::Approx(1e5));
  REQUIRE(con.ghost().flow_state()(0)(2)[4] == Catch::Approx(3.4));
}

TEST_CASE("Nonpenetration") {
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 4, 2, row_size};
  hexed::Tree tree(3, 1.);
  hexed::Deformed_element element {params, tree};
  hexed::Nonpenetration nonpen;
  hexed::Boundary_connection con(element.face(2), 0, 0);
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    double qpoint_nrml [] {-4., 3.};
    for (int i_dim = 0; i_dim < 2; ++i_dim) con.normal()(i_dim)[i_qpoint] = qpoint_nrml[i_dim];
    double state [] {1., 1., 1.2, 1e5/0.4 + 0.5*1.2*2.};
    for (int i = 0; i < 2; ++i) {
      for (int i_var = 0; i_var < 4; ++i_var) con.inside().flow_state()(i)(i_var)[i_qpoint] = state[i_var];
    }
  }
  SECTION("apply_state") {
    nonpen.apply_state(con);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      // require tangential momentum unchanged
      REQUIRE(  3*con.ghost().flow_state()(0)(0)[i_qpoint]
              + 4*con.ghost().flow_state()(0)(1)[i_qpoint] == Catch::Approx(7.));
      // require normal momentum flipped
      REQUIRE( -4*con.ghost().flow_state()(0)(0)[i_qpoint]
              + 3*con.ghost().flow_state()(0)(1)[i_qpoint] == Catch::Approx(1.));
    }
  }
  SECTION("apply_flux") {
    nonpen.apply_flux(con);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      // require tangential momentum flux flipped
      REQUIRE(  3*con.ghost().flow_state()(1)(0)[i_qpoint]
              + 4*con.ghost().flow_state()(1)(1)[i_qpoint] == Catch::Approx(-7.));
      // require normal momentum flux unchanged
      REQUIRE( -4*con.ghost().flow_state()(1)(0)[i_qpoint]
              + 3*con.ghost().flow_state()(1)(1)[i_qpoint] == Catch::Approx(-1.));
      // require scalar flux flipped
      REQUIRE(con.ghost().flow_state()(1)(2)[i_qpoint] == Catch::Approx(-con.inside().flow_state()(1)(2)[i_qpoint]));
      REQUIRE(con.ghost().flow_state()(1)(3)[i_qpoint] == Catch::Approx(-con.inside().flow_state()(1)(3)[i_qpoint]));
    }
  }
}

TEST_CASE("No_slip") {
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 4, 2, row_size};
  hexed::Tree tree(3, 1.);
  hexed::Deformed_element element {params, tree};
  hexed::Boundary_connection con(element.face(0), 0, 3);
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    for (int i_dim = 0; i_dim < 2; ++i_dim) con.normal()(i_dim)[i_qpoint] = .7/std::sqrt(2.);
  }
  double state [] {1., 1., 1.2, 1e5/0.4 + 0.5*1.2*2.};
  double flux [] {10., -20., 1.3, 10.};
  SECTION("isothermal") {
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) con.inside().flow_state()(0)(i_var)[i_qpoint] = state[i_var];
    }
    hexed::No_slip no_slip(std::make_shared<hexed::Prescribed_energy>(1e6), 1.4, hexed::inviscid, hexed::laminar);
    no_slip.apply_state(con);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_dim = 0; i_dim < 2; ++i_dim) {
        REQUIRE(con.ghost().flow_state()(0)(i_dim)[i_qpoint] == Catch::Approx(-1.));
      }
      REQUIRE(con.ghost().flow_state()(0)(2)[i_qpoint] == Catch::Approx(1.2));
      REQUIRE(std::sqrt(con.ghost().flow_state()(0)(3)[i_qpoint]*con.inside().flow_state()(0)(3)[i_qpoint])
              == Catch::Approx(1e6*1.2));
    }
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) con.inside().flow_state()(1)(i_var)[i_qpoint] = flux[i_var];
    }
    no_slip.apply_flux(con);
    flux[2] *= -1; // mass flux should always inverted
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(con.ghost().flow_state()(1)(i_var)[i_qpoint] == Catch::Approx(flux[i_var]));
      }
    }
  }
  SECTION("specified flux") {
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) con.inside().flow_state()(0)(i_var)[i_qpoint] = state[i_var];
    }
    hexed::No_slip no_slip(std::make_shared<hexed::Prescribed_heat_flux>(3.),
                           1.4, hexed::inviscid, hexed::laminar);
    no_slip.apply_state(con);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_dim = 0; i_dim < 2; ++i_dim) REQUIRE(con.ghost().flow_state()(0)(i_dim)[i_qpoint] == Catch::Approx(-1.));
      REQUIRE(con.ghost().flow_state()(0)(2)[i_qpoint] == Catch::Approx(state[2]));
      REQUIRE(con.ghost().flow_state()(0)(3)[i_qpoint] == Catch::Approx(state[3]));
    }
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) con.inside().flow_state()(1)(i_var)[i_qpoint] = flux[i_var];
    }
    no_slip.apply_flux(con);
    flux[2] *= -1; // mass flux should always inverted
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 3; ++i_var) {
        REQUIRE(con.ghost().flow_state()(1)(i_var)[i_qpoint] == Catch::Approx(flux[i_var]));
      }
      REQUIRE((con.ghost().flow_state()(1)(3)[i_qpoint] + con.inside().flow_state()(1)(3)[i_qpoint])/2 == Catch::Approx(-3.*.7));
    }
  }
  SECTION("specified emissivity") {
    state[3] = 1e5/.4;
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) {
        con.inside().flow_state()(0)(i_var)[i_qpoint] = state[i_var];
        con.inside().flow_state()(1)(i_var)[i_qpoint] = flux[i_var];
      }
    }
    auto thermal = std::make_shared<hexed::Thermal_equilibrium>();
    thermal->emissivity = .8;
    hexed::No_slip no_slip(thermal, 1.4, hexed::inviscid, hexed::laminar);
    double temp = 1e5/1.2/hexed::constants::specific_gas_air;
    no_slip.apply_state(con);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      for (int i_var = 0; i_var < 4; ++i_var) {
        con.inside().flow_state()(0)(i_var)[i_qpoint] = 0;
      }
    }
    no_slip.apply_flux(con);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
      REQUIRE((con.ghost().flow_state()(1)(3)[i_qpoint] + con.inside().flow_state()(1)(3)[i_qpoint])/2 == Catch::Approx(-.8*hexed::constants::stefan_boltzmann*std::pow(temp, 4)*.7));
    }
  }
}

TEST_CASE("Copy") {
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 5, 3, row_size};
  hexed::Tree tree(3, 1.);
  hexed::Element element {params, tree};
  const int n_qpoint = row_size*row_size;
  hexed::Copy copy;
  hexed::Boundary_connection con(element.face(2), 0, 0);
  // set inside face to something arbitrary
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    con.inside().flow_state()(0)(0)[i_qpoint] = 20.;
    con.inside().flow_state()(0)(1)[i_qpoint] = -10.;
    con.inside().flow_state()(0)(2)[i_qpoint] = 50.;
    con.inside().flow_state()(0)(3)[i_qpoint] = 0.9;
    con.inside().flow_state()(0)(4)[i_qpoint] = 1e4;
  }
  copy.apply_state(con);
  // check that ghost face is equal to inside
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    REQUIRE(con.ghost().flow_state()(0)(0)[i_qpoint] == Catch::Approx(20.));
    REQUIRE(con.ghost().flow_state()(0)(1)[i_qpoint] == Catch::Approx(-10.));
    REQUIRE(con.ghost().flow_state()(0)(2)[i_qpoint] == Catch::Approx(50.));
    REQUIRE(con.ghost().flow_state()(0)(3)[i_qpoint] == Catch::Approx(0.9));
    REQUIRE(con.ghost().flow_state()(0)(4)[i_qpoint] == Catch::Approx(1e4));
  }
}
