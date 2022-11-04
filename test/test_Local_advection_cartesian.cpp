#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Local_advection_cartesian.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Local_advection_cartesian")
{
  double d_t = 0.2;

  const int n_elem = 2;
  const int row_size = 2;
  hexed::Storage_params params {3, 4, 2, row_size};
  int n_qpoint = params.n_qpoint();
  std::vector<std::unique_ptr<hexed::Element>> elements;
  hexed::Gauss_legendre basis (row_size);
  /*
   * to test the correctness of the time derivative, set the RK weight to 0.5
   * and the RK reference state to minus the initial state. Thus the initial state
   * and the reference state cancel out, and the result is 0.5 times the time derivative
   * being written to the state.
   */
  double rk_weight = 0.5;

  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    elements.emplace_back(new hexed::Element {params, {}, 2.});
    double* state = elements[i_elem]->stage(0);
    double* tss = elements[i_elem]->time_step_scale();
    double* face = elements[i_elem]->face();

    #define SET_VARS \
      double veloc [] {2, pos1}; \
      double scalar = pos0; \

    for (int i = 0; i < row_size; ++i)
    {
      for (int j = 0; j < row_size; ++j)
      {
        int i_qpoint = i*row_size + j;
        double pos0 = 2*basis.node(i); double pos1 = 2*basis.node(j);
        SET_VARS
        // set initial state
        state[0*n_qpoint + i_qpoint] = veloc[0];
        state[1*n_qpoint + i_qpoint] = veloc[1];
        state[2*n_qpoint + i_qpoint] = scalar;
        // set RK reference state to negative of initial state
        state[3*n_qpoint + i_qpoint] = -state[2*n_qpoint + i_qpoint];
        tss[i_qpoint] = .7; // make sure time step scale has not effect
      }
      // set face state to match interior state
      for (int i_dim : {0, 1})
      {
        for (int positive : {0, 1})
        {
          double* qpoint_start = face + (2*i_dim + positive)*row_size*4 + i;
          double pos0 {2*(i_dim ? basis.node(i) : positive)};
          double pos1 {2*(i_dim ? positive : basis.node(i))};
          SET_VARS
          qpoint_start[2*row_size] = scalar*veloc[i_dim];
        }
      }
    }
    #undef SET_VARS
  }
  car_elem_view elem_view {elements};
  (*hexed::kernel_factory<hexed::Local_advection_cartesian>(2, row_size, basis, d_t, rk_weight))(elem_view);
  for (auto& element : elements) {
    double* state = element->stage(0);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      double pos0 = 2*basis.node(i_qpoint/row_size);
      REQUIRE(state[2*n_qpoint + i_qpoint] == Approx(-.5*.2*(2 + pos0)).scale(1.));
    }
  }
}