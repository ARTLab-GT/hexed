#include <catch2/catch.hpp>
#include <config.hpp>
#include <Surface_rotation.hpp>

TEST_CASE("Surface_rotation")
{
  const int n_dim = 2;
  const int row_size = std::min(7, hexed::config::max_row_size);
  const int n_qpoint = row_size;
  double state[(n_dim+2)*n_qpoint];
  double jacobian[n_dim*n_dim*n_qpoint];
  double jacobian_point1 [n_dim*n_dim] {3., 0., 0., 1.};
  double jacobian_rest [n_dim*n_dim] {2., 1., -2., 0.};
  const int point1 = 3;
  double state_all [n_dim+2] {2., 2., 1.3, 1e4};
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) {
      jacobian[i_jac*n_qpoint + i_qpoint] = jacobian_rest[i_jac];
    }
    for (int i_var = 0; i_var < n_dim+2; ++i_var) state[i_var*n_qpoint + i_qpoint] = state_all[i_var];
  }
  for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) {
    jacobian[i_jac*n_qpoint + point1] = jacobian_point1[i_jac];
  }
  auto surf_rot {hexed::kernel_factory<hexed::Surface_rotation>(n_dim, row_size, (double*)jacobian, 1)};

  // test that the transformed normal component is correct and Jacobian determinant are correct
  surf_rot->to_surface(state);
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    if (i_qpoint == point1) {
      REQUIRE(state[1*n_qpoint + i_qpoint] == Approx(2.));
      REQUIRE(surf_rot->jacobian_determinant(i_qpoint) == Approx(3.));
    } else {
      REQUIRE(state[1*n_qpoint + i_qpoint] == Approx(2.*std::sqrt(2.)));
      REQUIRE(surf_rot->jacobian_determinant(i_qpoint) == Approx(2.*std::sqrt(2.)));
    }
  }

  // test that the inverse transformation takes everything back to normal
  surf_rot->from_surface(state);
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    for (int i_var = 0; i_var < n_dim+2; ++i_var) {
      REQUIRE(state[i_var*n_qpoint + i_qpoint] == Approx(state_all[i_var]));
    }
  }
}
