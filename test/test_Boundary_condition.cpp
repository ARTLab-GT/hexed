#include <catch2/catch.hpp>
#include <config.hpp>
#include <Boundary_condition.hpp>

class Dummy : public hexed::Boundary_condition
{
  public:
  virtual void apply(hexed::Boundary_face&) {}
};

TEST_CASE("Typed_boundary_connection")
{
  hexed::Storage_params params {3, 4, 2, 4};
  hexed::Element element {params};
  Dummy bc;
  hexed::Typed_bound_connection<hexed::Element> tbc0 {element, 1, false, 0};
  hexed::Typed_bound_connection<hexed::Element> tbc1 {element, 1,  true, 1};
  REQUIRE(tbc0.storage_params().n_var == 4);
  // check that the correct face of the element is retrieved
  REQUIRE(tbc0.inside_face() == element.face() + (2*1 + 0)*4*4);
  REQUIRE(tbc1.inside_face() == element.face() + (2*1 + 1)*4*4);
  // check that ghost data exists (otherwise segfault)
  tbc0.ghost_face()[0] = 1.;
  tbc0.ghost_face()[4*4 - 1] = 1.;
  // check that order of faces is correct
  REQUIRE(tbc0.face(0) == tbc0.inside_face());
  REQUIRE(tbc0.face(1) == tbc0.ghost_face());
  REQUIRE(tbc1.face(0) == tbc1.inside_face());
  REQUIRE(tbc1.face(1) == tbc1.ghost_face());
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
  hexed::Freestream freestream {{10., 30., -20., 1.3, 1.2e5}};
  hexed::Typed_bound_connection<hexed::Element> tbc {element, 1, false, 0};
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
  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 4, 2, row_size};
  hexed::Deformed_element element {params};
  hexed::Nonpenetration nonpen;
  hexed::Typed_bound_connection<hexed::Deformed_element> tbc {element, 0, true, 0};
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    double qpoint_jacobian [] {1., 3., 0., 4.};
    for (int i_jac = 0; i_jac < 4; ++i_jac) tbc.jacobian()[i_jac*row_size + i_qpoint] = qpoint_jacobian[i_jac];
    double state [] {1., 1., 1.2, 1e5/0.4 + 0.5*1.2*2.};
    for (int i_var = 0; i_var < 4; ++i_var) tbc.inside_face()[i_var*row_size + i_qpoint] = state[i_var];
  }
  nonpen.apply(tbc);
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    // require tangential momentum unchanged
    REQUIRE(  3*tbc.ghost_face()[0*row_size + i_qpoint]
            + 4*tbc.ghost_face()[1*row_size + i_qpoint] == Approx(7.));
    REQUIRE( -4*tbc.ghost_face()[0*row_size + i_qpoint]
            + 3*tbc.ghost_face()[1*row_size + i_qpoint] == Approx(1.));
    // require normal momentum flipped
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
  copy.apply(tbc);
  // check that ghost face is equal to inside
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    REQUIRE(tbc.ghost_face()[0*n_qpoint + i_qpoint] == Approx(20.));
    REQUIRE(tbc.ghost_face()[1*n_qpoint + i_qpoint] == Approx(-10.));
    REQUIRE(tbc.ghost_face()[2*n_qpoint + i_qpoint] == Approx(50.));
    REQUIRE(tbc.ghost_face()[3*n_qpoint + i_qpoint] == Approx(0.9));
    REQUIRE(tbc.ghost_face()[4*n_qpoint + i_qpoint] == Approx(1e4));
  }
}
