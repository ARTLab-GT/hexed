#include <catch2/catch.hpp>

#include <Deformed_element.hpp>

void assert_equal(std::array<double, 3> computed, std::array<double, 3> correct)
{
  for (int i_dim = 0; i_dim < 3; ++i_dim)
  {
    CHECK(computed[i_dim] == Approx(correct[i_dim]).margin(1e-14));
  }
}

TEST_CASE("Deformed_element.hpp")
{
  cartdg::Storage_params params {2, 2, 2, 4};
  cartdg::Deformed_element element {params};
  for (int i_jac = 0; i_jac < 2*2*16; ++i_jac)
  {
    element.jacobian()[i_jac] = i_jac/4.;
  }
  for (int i_dim = 0, i_jac = 0; i_dim < 2; ++i_dim)
  {
    for (int j_dim = 0; j_dim < 2; ++j_dim)
    {
      for (int i_qpoint = 0; i_qpoint < 16; ++i_qpoint)
      {
        REQUIRE(element.jacobian(i_dim, j_dim, i_qpoint) == i_jac++/4.);
      }
    }
  }

  element.jacobian()[0*16 + 15] = 1.;
  element.jacobian()[1*16 + 15] = 2.;
  element.jacobian()[2*16 + 15] = .1;
  element.jacobian()[3*16 + 15] = 3.;
  REQUIRE(element.jacobian_determinant(15) == Approx(2.8));

  for (int i_adj = 0; i_adj < 16; ++i_adj)
  {
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
  for (int i_jac = 0; i_jac < 9; ++i_jac)
  {
    element3d.jacobian()[i_jac*8] = jacobian[i_jac];
  }
  REQUIRE(element3d.jacobian_determinant(0) == Approx(1.*.7 + 0.3*1.5));

  SECTION("vertex arrangement")
  {
    // check that vertices start out in correct location
    cartdg::Storage_params params3d {3, 5, 3, 4};
    cartdg::Deformed_element elem3d {params3d, {1, 2, -1}, 0.2};
    assert_equal(elem3d.vertex(0).pos, {0.2, 0.4, -0.2});
    assert_equal(elem3d.vertex(1).pos, {0.2, 0.4,  0. });
    assert_equal(elem3d.vertex(2).pos, {0.2, 0.6, -0.2});
    assert_equal(elem3d.vertex(5).pos, {0.4, 0.4,  0. });
    assert_equal(elem3d.vertex(7).pos, {0.4, 0.6,  0. });

    // check that vertices have correct neighbors
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
