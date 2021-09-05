#include <catch.hpp>

#include <Element.hpp>

TEST_CASE("Element")
{
  cartdg::Storage_params params {4, 5, 3, 6};
  unsigned n_dof = params.n_dof();
  cartdg::Element element {params};
  for (unsigned i_stage = 0; i_stage < 4; ++i_stage)
  {
    for (unsigned i_dof = 0; i_dof < n_dof; ++i_dof)
    {
      element.stage(i_stage)[i_dof] = 0.;
    }
  }
  for (unsigned i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(0)[i_dof] = 1.2;
  for (unsigned i_dof = 0; i_dof < n_dof; ++i_dof) REQUIRE(element.stage(3)[i_dof] == 0.);
  for (unsigned i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(3)[i_dof] = 1.3;
  for (unsigned i_dof = 0; i_dof < n_dof; ++i_dof)
  {
    REQUIRE(element.stage(0)[i_dof] == 1.2);
    REQUIRE(element.stage(1)[i_dof] == 0.);
    REQUIRE(element.stage(2)[i_dof] == 0.);
    REQUIRE(element.stage(3)[i_dof] == 1.3);
  }

  cartdg::Element copy {element};
  cartdg::Element assigned {params};
  assigned = element;
  element.stage(1)[0] = 1.4;
  for (cartdg::Element* test_elem : {&copy, &assigned})
  {
    for (unsigned i_dof = 0; i_dof < n_dof; ++i_dof)
    {
      REQUIRE(test_elem->stage(0)[i_dof] == 1.2);
      REQUIRE(test_elem->stage(1)[i_dof] == 0.);
      REQUIRE(test_elem->stage(2)[i_dof] == 0.);
      REQUIRE(test_elem->stage(3)[i_dof] == 1.3);
    }
  }
}
