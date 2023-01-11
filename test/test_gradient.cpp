#include <catch2/catch.hpp>
#include <hexed/Spatial.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/pde.hpp>

TEST_CASE("Spatial::gradient")
{
  constexpr int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {1, 4, 2, row_size};
  hexed::Deformed_element elem(params, {}, .8);
  elem.vertex(3).pos = {.9, .7, 0.};
  double face_data[4][4*row_size];
  for (int i_face = 0; i_face < 4; ++i_face) elem.faces[i_face] = face_data[i_face];
  hexed::Gauss_legendre basis(row_size);
  elem.set_jacobian(basis);
  hexed::Mat<2, 3> coefs;
  coefs << .1, -.4, .3,
           -.2, .6, .6;
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
    for (int i_var = 0; i_var < 3; ++i_var) {
      auto pos = elem.position(basis, i_qpoint);
      double total = 0;
      for (int i_dim = 0; i_dim < 2; ++i_dim) total += coefs(i_dim, i_var)*pos[i_dim];
      elem.stage(0)[i_var*params.n_qpoint() + i_qpoint] = std::exp(total);
    }
  }
  double gradient[2][3][row_size*row_size];
  hexed::Derivative<row_size> deriv(basis);
  hexed::Spatial<hexed::Deformed_element, hexed::pde::Euler, true>::compute_gradient<2, row_size, 3>(elem.stage(0), elem.reference_level_normals(), gradient[0][0], deriv);
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
    for (int i_var = 0; i_var < 3; ++i_var) {
      auto pos = elem.position(basis, i_qpoint);
      double total = 0;
      for (int i_dim = 0; i_dim < 2; ++i_dim) total += coefs(i_dim, i_var)*pos[i_var];
      for (int i_dim = 0; i_dim < 2; ++i_dim) {
        REQUIRE(gradient[i_dim][i_var][i_qpoint] == Approx(coefs(i_dim, i_var)*elem.stage(0)[i_var*params.n_qpoint() + i_qpoint]).epsilon(1e-4));
      }
    }
  }
}
