#include <catch2/catch.hpp>
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
  cartdg::Typed_bound_connection<cartdg::Element> tbc0 {element, {1}, bc, false};
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
