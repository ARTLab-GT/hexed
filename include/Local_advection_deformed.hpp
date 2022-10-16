#ifndef HEXED_LOCAL_ADVECTION_DEFORMED_HPP_
#define HEXED_LOCAL_ADVECTION_DEFORMED_HPP_

#include "Vector_view.hpp"
#include "Deformed_element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Write_face.hpp"

namespace hexed
{

/*
 * Computes the local update for one Runge-Kutta stage of linear advection.
 * See `Local_advection_cartesian`.
 */
template <int n_dim, int row_size>
class Local_advection_deformed : public Kernel<Deformed_element&>
{
  Derivative<row_size> derivative;
  Write_face<n_dim, row_size> write_face;
  double dt;
  double rkw;

  public:
  Local_advection_deformed(const Basis& basis, double d_time, double rk_weight) :
    derivative{basis},
    write_face{basis},
    dt{d_time},
    rkw{rk_weight}
  {}

  virtual void operator()(Sequence<Deformed_element&>& elements)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    constexpr int n_face_dof = n_var*n_qpoint/row_size;

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      auto& elem {elements[i_elem]};
      double* veloc = elements[i_elem].stage(0); // advection velocity, not physical velocity
      double* state = veloc + n_dim*n_qpoint;
      double* rk_reference = state + n_var*n_qpoint;
      double time_rate [n_qpoint] {};
      double* normals = elem.reference_level_normals();
      double* determinant = elem.jacobian_determinant();
      double* face = elem.face();
      const double d_t_by_d_pos = dt/elem.nominal_size();

      // Compute update
      for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
           stride /= row_size, n_rows *= row_size, ++i_dim)
      {
        double* dim_nrml = normals + i_dim*n_dim*n_qpoint;
        int i_face_qpoint {0};
        for (int i_outer = 0; i_outer < n_rows; ++i_outer)
        {
          for (int i_inner = 0; i_inner < stride; ++i_inner)
          {
            // compute velocity in reference space by taking dot product with normals
            double row_v [row_size] {};
            for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                int i = j_dim*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride;
                row_v[i_qpoint] += veloc[i]*dim_nrml[i];
              }
            }
            // fetch state
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
        double updated = time_rate[i_qpoint]/determinant[i_qpoint]*d_t_by_d_pos + state[i_qpoint];
        state[i_qpoint] = rkw*updated + (1. - rkw)*rk_reference[i_qpoint];
      }
      write_face(state, face);
    }
  }
};

template<>
class Kernel_traits<Local_advection_deformed>
{
  public:
  using base_t = Kernel<Deformed_element&>;
};

}
#endif
