#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <math.hpp>
#include <Deformed_element.hpp>
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
    cartdg::Element elem3d {params3d, {1, 2, -1}, 0.2};
    elem3d.vertex(5).mobile = true;
    elem3d.vertex(5).calc_relax();
    elem3d.vertex(5).apply_relax();
    assert_equal(elem3d.vertex(5).pos, {0.4 - 0.1/3., 0.4 + 0.1/3., -0.1/3.});

    cartdg::Storage_params params2d {3, 4, 2, 4};
    cartdg::Deformed_element elem2d {params2d, {1, 2, -1}, 1.};
    assert_equal(elem2d.vertex(0).pos, {1., 2., 0.});
    assert_equal(elem2d.vertex(1).pos, {1., 3., 0.});
    assert_equal(elem2d.vertex(2).pos, {2., 2., 0.});
    elem2d.vertex(0).mobile = true;
    elem2d.vertex(0).calc_relax();
    elem2d.vertex(0).apply_relax();
    assert_equal(elem2d.vertex(0).pos, {1.25, 2.25, 0.});

    cartdg::Deformed_element elem3d_1 {params3d, {3,}, 1.};
    assert_equal(elem3d_1.vertex(0).pos, {3., 0., 0.});
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
