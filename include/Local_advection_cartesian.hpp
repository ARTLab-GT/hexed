#ifndef HEXED_LOCAL_ADVECTION_CARTESIAN_HPP_
#define HEXED_LOCAL_ADVECTION_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Write_face.hpp"

namespace hexed
{

/*
 * Computes the local update for the advection equation for one Runge-Kutta stage.
 * For the face data and stage 0 of the interior data, the first `n_dim` variables
 * should contain the advection velocity
 * and variable `n_dim` should contain the advected scalar.
 * Variable `n_dim + 1` of stage 0 should contain the RK reference for the scalar.
 */
template <int n_dim, int row_size>
class Local_advection_cartesian : public Kernel<Element&>
{
  Derivative<row_size> derivative;
  Write_face<n_dim, row_size> write_face;
  double update;
  double curr;
  double ref;

  public:
  Local_advection_cartesian(const Basis& basis, double update_coef, double current_coef, double reference_coef) :
    derivative{basis},
    write_face{basis},
    update{update_coef},
    curr{current_coef},
    ref{reference_coef}
  {}

  virtual void operator()(Sequence<Element&>& elements)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    constexpr int n_face_dof = n_var*n_qpoint/row_size;

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      double* veloc = elements[i_elem].stage(0); // advection velocity, not physical velocity
      double* state = veloc + n_dim*n_qpoint;
      double* rk_reference = state + n_qpoint;
      double time_rate [n_qpoint] {};
      double* face = elements[i_elem].face();
      const double d_pos = elements[i_elem].nominal_size();

      // Compute update
      for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
           stride /= row_size, n_rows *= row_size, ++i_dim)
      {
        int i_face_qpoint {0};
        for (int i_outer = 0; i_outer < n_rows; ++i_outer)
        {
          for (int i_inner = 0; i_inner < stride; ++i_inner)
          {
            // Fetch this row of data
            double row_v [row_size];
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              row_v[i_qpoint] = veloc[i_dim*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
            double row_s [row_size];
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              row_s[i_qpoint] = state[i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
            // Calculate flux
            Eigen::Matrix<double, row_size, 1> flux;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              flux(i_qpoint) = row_v[i_qpoint]*row_s[i_qpoint];
            }
            // fetch boundary flux
            Eigen::Matrix<double, 2, 1> boundary_values;
            for (int is_positive : {0, 1}) {
              boundary_values(is_positive) = face[(i_dim*2 + is_positive)*n_face_dof + n_dim*n_qpoint/row_size + i_face_qpoint];
            }
            // Differentiate flux
            Eigen::Matrix<double, row_size, 1> row_w = -derivative(flux, boundary_values);
            // Add dimensional component to update
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              int offset = i_outer*stride*row_size + i_inner + i_qpoint*stride;
              time_rate[offset] += row_w(i_qpoint);
            }
            ++i_face_qpoint;
          }
        }
      }

      // write the updated solution
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        state[i_qpoint] = update*time_rate[i_qpoint]/d_pos
                          + curr*state[i_qpoint]
                          + ref*rk_reference[i_qpoint];
      }
      write_face(veloc, face);
    }
  }
};

}
#endif
