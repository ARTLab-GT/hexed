#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Local_advection_deformed.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Local_advection_deformed")
{
  double d_t = 0.2;

  const int row_size = hexed::config::max_row_size;
  hexed::Storage_params params {3, 4, 2, row_size};
  int n_qpoint = params.n_qpoint();
  std::vector<std::unique_ptr<hexed::Deformed_element>> elements;
  hexed::Gauss_legendre basis (row_size);
  /*
   * to test the correctness of the time derivative, set the RK weight to 0.5
   * and the RK reference state to minus the initial state. Thus the initial state
   * and the reference state cancel out, and the result is 0.5 times the time derivative
   * being written to the state.
   */
  double rk_weight = 0.5;

  elements.emplace_back(new hexed::Deformed_element {params, {}, 2.});
  auto& elem = *elements.back();
  // deform one of the vertices
  elem.vertex(1).pos = {.05*2., -.1*2., 0.};
  elem.set_jacobian(basis);
  double* state = elem.stage(0);
  double* face  = elem.face();
  double face_normal [4][2*n_qpoint];
  for (int i_face = 0; i_face < 4; ++i_face) {
    elem.face_normal(i_face) = face_normal[i_face];
    for (int i_data = 0; i_data < n_qpoint/row_size*2; ++i_data) {
      face_normal[i_face][i_data] = elem.face()[i_face*n_qpoint/row_size*4 + i_data];
    }
  }

  #define SET_VARS \
    double veloc [] {2, pos[1]}; \
    double scalar = pos[0]; \

  for (int i = 0; i < row_size; ++i)
  {
    for (int j = 0; j < row_size; ++j)
    {
      int i_qpoint = i*row_size + j;
      auto pos = elem.position(basis, i_qpoint);
      SET_VARS
      // set initial state
      state[0*n_qpoint + i_qpoint] = veloc[0];
      state[1*n_qpoint + i_qpoint] = veloc[1];
      state[2*n_qpoint + i_qpoint] = scalar;
      // set RK reference state to negative of initial state
      state[3*n_qpoint + i_qpoint] = -state[2*n_qpoint + i_qpoint];
    }
  }
  #undef SET_VARS
  {
    hexed::Vector_view<hexed::Element&, std::unique_ptr<hexed::Deformed_element>, &hexed::ptr_convert<hexed::Element&, std::unique_ptr<hexed::Deformed_element>>> elem_view {elements};
    (*hexed::kernel_factory<hexed::Write_face>(2, row_size, basis))(elem_view);
  }
  // set face state to match interior state
  for (int i = 0; i < row_size; ++i) {
    for (int i_dim : {0, 1}) {
      for (int positive : {0, 1}) {
        double flux = 0;
        double* qpoint_start = face + (2*i_dim + positive)*row_size*4 + i;
        for (int j_dim = 0; j_dim < 2; ++j_dim) {
          flux += face_normal[2*i_dim + positive][j_dim*row_size + i]
                  *qpoint_start[j_dim*row_size]*qpoint_start[2*row_size];
        }
        qpoint_start[2*row_size] = flux;
      }
    }
  }
  {
    def_elem_view elem_view {elements};
    (*hexed::kernel_factory<hexed::Local_advection_deformed>(2, row_size, basis, d_t, rk_weight))(elem_view);
  }
  for (auto& element : elements) {
    double* state = element->stage(0);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      double pos0 = elem.position(basis, i_qpoint)[0];
      REQUIRE(state[2*n_qpoint + i_qpoint] == Approx(-.5*.2*(2 + pos0)).scale(1.));
    }
  }
}
