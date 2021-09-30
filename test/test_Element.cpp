#include <catch2/catch.hpp>

#include <Element.hpp>

TEST_CASE("Element")
{
  cartdg::Storage_params params {4, 5, 3, 6};
  int n_dof = params.n_dof();
  cartdg::Element element {params};
  // sometimes we will write to and read from some data just
  // to be sure the storage is there and doesn't overlap
  for (int i_stage = 0; i_stage < 4; ++i_stage)
  {
    for (int i_dof = 0; i_dof < n_dof; ++i_dof)
    {
      element.stage(i_stage)[i_dof] = 0.;
    }
  }
  for (int i_face_data = 0; i_face_data < n_dof/params.row_size*3*2; ++i_face_data)
  {
    element.face()[i_face_data] = 1.;
  }
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(0)[i_dof] = 1.2;
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) REQUIRE(element.stage(3)[i_dof] == 0.);
  for (int i_dof = 0; i_dof < n_dof; ++i_dof) element.stage(3)[i_dof] = 1.3;
  for (int i_dof = 0; i_dof < n_dof; ++i_dof)
  {
    REQUIRE(element.stage(0)[i_dof] == 1.2);
    REQUIRE(element.stage(1)[i_dof] == 0.);
    REQUIRE(element.stage(2)[i_dof] == 0.);
    REQUIRE(element.stage(3)[i_dof] == 1.3);
  }
  for (int i_face_data = 0; i_face_data < n_dof/params.row_size*3*2; ++i_face_data)
  {
    REQUIRE(element.face()[i_face_data] == 1.);
  }
  for (int i_vert = 0; i_vert < 8; ++i_vert)
  {
    REQUIRE(element.viscosity()[i_vert] == 0.);
  }
  REQUIRE(element.viscous() == false);
  element.viscosity()[3] = 0.1;
  REQUIRE(element.viscous() == true);
  element.derivative()[0] = 1.;
  element.derivative()[6*6*6 - 1] = 1.;
  REQUIRE(element.derivative()[0] == 1.);
  REQUIRE(element.derivative()[6*6*6 - 1] == 1.);
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint)
  {
    REQUIRE(element.jacobian(0, 0, i_qpoint) == 1.);
    REQUIRE(element.jacobian(1, 0, i_qpoint) == 0.);
    REQUIRE(element.jacobian(0, 2, i_qpoint) == 0.);
    REQUIRE(element.jacobian(2, 2, i_qpoint) == 1.);
    REQUIRE(element.jacobian_determinant(i_qpoint) == 1.);
  }

  cartdg::Element copy {element};
  cartdg::Element assigned {params};
  assigned = element;
  element.stage(1)[0] = 1.4;
  for (cartdg::Element* test_elem : {&copy, &assigned})
  {
    for (int i_dof = 0; i_dof < n_dof; ++i_dof)
    {
      REQUIRE(test_elem->stage(0)[i_dof] == 1.2);
      REQUIRE(test_elem->stage(1)[i_dof] == 0.);
      REQUIRE(test_elem->stage(2)[i_dof] == 0.);
      REQUIRE(test_elem->stage(3)[i_dof] == 1.3);
    }
  }
}
