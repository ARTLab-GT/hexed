#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <math.hpp>
#include <Deformed_element.hpp>
#include <Equidistant.hpp>
#include <Gauss_legendre.hpp>
#include "testing_utils.hpp"

TEST_CASE("Deformed_element")
{
  cartdg::Storage_params params {2, 2, 2, 4};
  cartdg::Deformed_element element {params};
  for (int i_jac = 0; i_jac < 2*2*16; ++i_jac) {
    element.jacobian()[i_jac] = i_jac/4.;
  }
  for (int i_dim = 0, i_jac = 0; i_dim < 2; ++i_dim) {
    for (int j_dim = 0; j_dim < 2; ++j_dim) {
      for (int i_qpoint = 0; i_qpoint < 16; ++i_qpoint) {
        REQUIRE(element.jacobian(i_dim, j_dim, i_qpoint) == i_jac++/4.);
      }
    }
  }

  element.jacobian()[0*16 + 15] = 1.;
  element.jacobian()[1*16 + 15] = 2.;
  element.jacobian()[2*16 + 15] = .1;
  element.jacobian()[3*16 + 15] = 3.;
  REQUIRE(element.jacobian_determinant(15) == Approx(2.8));

  for (int i_adj = 0; i_adj < 16; ++i_adj) {
    REQUIRE(element.node_adjustments()[i_adj] == 0.);
  }
  element.node_adjustments()[0] = 2.7;
  element.node_adjustments()[4*2*2 - 1] = 0.04;
  REQUIRE(element.node_adjustments()[0] == 2.7);
  REQUIRE(element.node_adjustments()[4*2*2 - 1] == 0.04);
  REQUIRE(element.vertex(1).mobile);

  cartdg::Storage_params params3d {1, 1, 3, 2};
  cartdg::Deformed_element element3d {params3d};
  double jacobian [9] {2., 1., -.7,
                       .5, 1.,  0.,
                       1., 0.,  .3};
  for (int i_jac = 0; i_jac < 9; ++i_jac) {
    element3d.jacobian()[i_jac*8] = jacobian[i_jac];
  }
  REQUIRE(element3d.jacobian_determinant(0) == Approx(1.*.7 + 0.3*1.5));

  SECTION("vertex relaxation")
  {
    cartdg::Storage_params params3d {3, 5, 3, 4};
    cartdg::Deformed_element elem3d {params3d, {1, 2, -1}, 0.2};
    elem3d.vertex(5).calc_relax();
    elem3d.vertex(5).apply_relax();
    assert_equal(elem3d.vertex(5).pos, {0.4 - 0.1/3., 0.4 + 0.1/3., -0.1/3.});

    cartdg::Storage_params params2d {3, 4, 2, 4};
    cartdg::Deformed_element elem2d {params2d, {1, 2, -1}, 1.};
    assert_equal(elem2d.vertex(0).pos, {1., 2., 0.});
    assert_equal(elem2d.vertex(1).pos, {1., 3., 0.});
    assert_equal(elem2d.vertex(2).pos, {2., 2., 0.});
    elem2d.vertex(0).calc_relax();
    elem2d.vertex(0).apply_relax();
    assert_equal(elem2d.vertex(0).pos, {1.25, 2.25, 0.});

    cartdg::Deformed_element elem3d_1 {params3d, {3,}, 1.};
    assert_equal(elem3d_1.vertex(0).pos, {3., 0., 0.});
  }

  SECTION("position calculation")
  {
    const int row_size = 3;
    static_assert (row_size <= cartdg::config::max_row_size);
    cartdg::Equidistant basis {row_size};
    cartdg::Storage_params params2 {2, 4, 2, row_size};
    cartdg::Deformed_element elem {params2};
    elem.vertex(3).pos[0] = 0.6;
    elem.node_adjustments()[2*3 + 1] =  0.2;
    elem.node_adjustments()[3*3 + 1] = -0.1;
    REQUIRE(elem.position(basis, 0)[0] == Approx(0.0));
    REQUIRE(elem.position(basis, 7)[0] == Approx(0.8));
    REQUIRE(elem.position(basis, 8)[0] == Approx(0.6));
    REQUIRE(elem.position(basis, 7)[1] == Approx(0.5));

    REQUIRE(elem.position(basis, 3)[0] == Approx(0.5 - 0.2*0.2));
    REQUIRE(elem.position(basis, 4)[0] == Approx(0.4 - 0.2*(0.2 - 0.1)/2));
    REQUIRE(elem.position(basis, 5)[0] == Approx(0.3 + 0.2*0.1));
    REQUIRE(elem.position(basis, 3)[1] == Approx(0.0 + 0.2));
    REQUIRE(elem.position(basis, 4)[1] == Approx(0.5 + (0.2 - 0.1)/2));
    REQUIRE(elem.position(basis, 5)[1] == Approx(1.0 - 0.1));
    // check that the face quadrature points are the same as the interior quadrature points
    // that happen to lie on the faces (true for equidistant and Lobatto bases but not Legendre)
    REQUIRE(elem.face_position(basis, 0, 2)[1] == elem.position(basis, 2)[1]);
    REQUIRE(elem.face_position(basis, 2, 1)[0] == elem.position(basis, 3)[0]);
    REQUIRE(elem.face_position(basis, 3, 1)[1] == elem.position(basis, 5)[1]);

    cartdg::Deformed_element elem1 {params2};
    elem1.node_adjustments()[1] = 0.1;
    REQUIRE(elem1.position(basis, 0)[0] == Approx(0.0));
    REQUIRE(elem1.position(basis, 6)[0] == Approx(1.0));
    REQUIRE(elem1.position(basis, 4)[0] == Approx(0.55));

    cartdg::Storage_params params3 {2, 5, 3, row_size};
    cartdg::Deformed_element elem2 {params3, {}, 0.2};
    elem2.node_adjustments()[4] = 0.01;
    REQUIRE(elem2.position(basis, 13)[0] == Approx(0.101));
    REQUIRE(elem2.position(basis, 13)[1] == Approx(.1));
    REQUIRE(elem2.position(basis, 13)[2] == Approx(.1));

    cartdg::Gauss_legendre leg_basis {row_size};
    cartdg::Deformed_element elem3 {params2, {}, 0.2};
    elem3.node_adjustments()[1] = 0.1;
    elem3.node_adjustments()[3] = -0.2;
    REQUIRE(elem3.position(leg_basis, 3)[0] == Approx(0.08));
    REQUIRE(elem3.position(leg_basis, 4)[0] == Approx(0.11));
  }

  SECTION("jacobian calculation")
  {
    const int row_size = 3;
    static_assert (row_size <= cartdg::config::max_row_size);
    cartdg::Equidistant basis {row_size};
    double* jac;
    cartdg::Storage_params params2 {2, 4, 2, row_size};
    cartdg::Deformed_element elem0 {params2, {0, 0}, 0.2};
    cartdg::Deformed_element elem1 {params2, {1, 1}, 0.2};
    elem0.vertex(3).pos = {0.8*0.2, 0.8*0.2, 0.};
    elem1.node_adjustments()[6 + 1] = 0.1;
    // jacobian is correct
    elem0.set_jacobian(basis);
    jac = elem0.jacobian();
    REQUIRE(jac[0*9    ] == Approx(1.));
    REQUIRE(jac[1*9    ] == Approx(0.));
    REQUIRE(jac[2*9    ] == Approx(0.));
    REQUIRE(jac[3*9    ] == Approx(1.));
    REQUIRE(jac[0*9 + 6] == Approx(1.));
    REQUIRE(jac[1*9 + 6] == Approx(-0.2));
    REQUIRE(jac[2*9 + 6] == Approx(0.));
    REQUIRE(jac[3*9 + 6] == Approx(0.8));
    REQUIRE(jac[0*9 + 8] == Approx(0.8));
    REQUIRE(jac[1*9 + 8] == Approx(-0.2));
    REQUIRE(jac[2*9 + 8] == Approx(-0.2));
    REQUIRE(jac[3*9 + 8] == Approx(0.8));
    elem1.set_jacobian(basis);
    jac = elem1.jacobian();
    REQUIRE(jac[0*9 + 5] == Approx(1.));
    REQUIRE(jac[1*9 + 5] == Approx(0.));
    REQUIRE(jac[2*9 + 5] == Approx(0.));
    REQUIRE(jac[3*9 + 5] == Approx(0.9));
    // time step is correct
    // At corner 3, diagonal is locally scaled by 1 - 2*(1 - 0.8) = 0.6 = (min singular value)
    REQUIRE(elem0.vertex_time_step_scale()[0] == 1.);
    REQUIRE(elem0.vertex_time_step_scale()[3] == Approx(0.6));

    cartdg::Storage_params params3 {2, 5, 3, row_size};
    cartdg::Deformed_element elem2 {params3, {0, 0, 0}, 0.2};
    elem2.vertex(7).pos = {0.8*0.2, 0.8*0.2, 0.8*0.2};
    elem2.set_jacobian(basis);
    jac = elem2.jacobian();
    REQUIRE(jac[0] == 1.);
    REQUIRE(jac[0*27 + 26] == Approx( 0.8));
    REQUIRE(jac[1*27 + 26] == Approx(-0.2));
    REQUIRE(jac[2*27 + 26] == Approx(-0.2));
    REQUIRE(jac[7*27 + 26] == Approx(-0.2));
    REQUIRE(jac[8*27 + 26] == Approx( 0.8));
  }
}

TEST_CASE("Deformed_face")
{
  const int n_dim = 3;
  const int row_size = std::min(4, cartdg::config::max_row_size);
  cartdg::Storage_params params {6, 7, n_dim, row_size}; // some more or less arbitrary parameters
  cartdg::Deformed_face face {params};
  face.jacobian()[0] = 0.8;
  face.jacobian()[n_dim*n_dim*cartdg::custom_math::pow(row_size, n_dim - 1) - 1] = 0.9;
  REQUIRE(face.jacobian(0, 0, 0) == 0.8);
  REQUIRE(face.jacobian(n_dim - 1, n_dim - 1, cartdg::custom_math::pow(row_size, n_dim - 1) - 1) == 0.9);
}

TEST_CASE("Deformed_elem_con")
{
  cartdg::Storage_params params {1, 3, 1, 1};
  cartdg::Deformed_element elem0 {params};
  cartdg::Deformed_element elem1 {params};
  cartdg::Face_index fi0 {&elem0, 2, 0};
  cartdg::Face_index fi1 {&elem1, 0, 1};
  cartdg::Deformed_elem_con con ({fi0, fi1});
  REQUIRE(con.face_index(0).element == &elem0);
  REQUIRE(con.face_index(0).i_dim == 2);
  REQUIRE(con.face_index(1).is_positive == 1);
  REQUIRE(con.flip_normal(0));
  REQUIRE(con.flip_normal(1));
  fi0.is_positive = 1;
  fi1.is_positive = 0;
  cartdg::Deformed_elem_con con1 ({fi0, fi1});
  REQUIRE(!con1.flip_normal(0));
  REQUIRE(!con1.flip_normal(1));
  SECTION("flip_tangential")
  {
    fi0.i_dim = 0;
    fi1.i_dim = 0;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(!con2.flip_tangential());
    }
    fi0.i_dim = 0;
    fi1.i_dim = 1;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(con2.flip_tangential());
    }
    fi0.i_dim = 1;
    fi1.i_dim = 0;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(con2.flip_tangential());
    }
    fi0.i_dim = 2;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(con2.flip_tangential());
    }
    fi1.is_positive = 1;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(!con2.flip_tangential());
    }
  }
  SECTION("transpose")
  {
    fi0.i_dim = 1;
    fi1.i_dim = 1;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(!con2.transpose());
    }
    fi1.i_dim = 0;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(!con2.transpose());
    }
    fi0.i_dim = 0;
    fi1.i_dim = 1;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(!con2.transpose());
    }
    fi0.i_dim = 1;
    fi1.i_dim = 2;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(!con2.transpose());
    }
    fi0.i_dim = 0;
    fi1.i_dim = 2;
    {
      cartdg::Deformed_elem_con con2 ({fi0, fi1});
      REQUIRE(con2.transpose());
    }
  }

  SECTION("position calculation")
  {
    const int row_size = 5;
    cartdg::Storage_params params {2, 4, 2, row_size};
    cartdg::Deformed_element elem {params};
    elem.vertex(3).pos[0] = 0.8;
    elem.vertex(3).pos[1] = 0.7;
    cartdg::Equidistant basis {row_size};
    // test first and last qpoints
    REQUIRE(elem.position(basis, 0).size() == 2);
    REQUIRE(elem.position(basis, 0)[0] == Approx(0.).scale(1.));
    REQUIRE(elem.position(basis, 0)[1] == Approx(0.).scale(1.));
    REQUIRE(elem.position(basis, params.n_qpoint() - 1)[0] == Approx(.8).scale(1.));
    REQUIRE(elem.position(basis, params.n_qpoint() - 1)[1] == Approx(.7).scale(1.));
    static_assert (row_size%2 == 1); // `row_size` must be odd for the following tests to work
    // test the qpoint at the midpoint of the positive-dimension0 face (the right-hand face)
    REQUIRE(elem.position(basis, row_size*(row_size - 1) + row_size/2)[0] == Approx(1. - 0.2/2).scale(1.));
    REQUIRE(elem.position(basis, row_size*(row_size - 1) + row_size/2)[1] == Approx(0.7/2).scale(1.));
    // test the qpoint at the middle of the element (the mean of the vertex positions)
    REQUIRE(elem.position(basis, params.n_qpoint()/2)[0] == Approx(.5 - 0.2/4).scale(1.));
    REQUIRE(elem.position(basis, params.n_qpoint()/2)[1] == Approx(.5 - 0.3/4).scale(1.));
  }
}

TEST_CASE("Deformed_elem_wall")
{
  cartdg::Storage_params params {1, 3, 1, 1};
  cartdg::Deformed_element elem {params};
  cartdg::Deformed_elem_wall dew {{&elem, 2, 1}, 4};
  REQUIRE(dew.face_index().element == &elem);
  REQUIRE(dew.face_index().i_dim == 2);
  REQUIRE(dew.face_index().is_positive);
  REQUIRE(dew.i_elem() == 4);
}
