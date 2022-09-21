#ifndef HEXED_LOCAL_DEFORMED_HPP_
#define HEXED_LOCAL_DEFORMED_HPP_

#include "Vector_view.hpp"
#include "Deformed_element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"

namespace hexed
{

/*
 * Computes the local update for one Runge-Kutta stage.
 * See `Local_cartesian_dynamic`.
 */
template <int n_dim, int row_size>
class Local_deformed : public Kernel<Deformed_element&>
{
  Derivative<row_size> derivative;
  double dt;
  double rkw;
  const double heat_rat;

  public:
  Local_deformed(const Basis& basis, double d_time, double rk_weight, double heat_ratio=1.4) :
    derivative{basis},
    dt{d_time},
    rkw{rk_weight},
    heat_rat{heat_ratio}
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
      double* state  = elem.stage(0);
      double* rk_reference = state + n_var*n_qpoint;
      double time_rate [n_var][n_qpoint] {};
      double* normals = elem.jacobian_data();
      double* determinant = normals + n_dim*n_dim*n_qpoint;
      double* face = elem.face();
      double* tss = elem.time_step_scale();
      const double d_t_by_d_pos = dt/elem.nominal_size();

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
            double row_r [n_var][row_size];
            double row_n [n_dim][row_size];
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                row_r[i_var][i_qpoint] = state[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
            for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                row_n[j_dim][i_qpoint] = normals[(i_dim*n_dim + j_dim)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
            // Calculate flux
            Eigen::Matrix<double, row_size, n_var> flux;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              #define READ(i) row_r[i][i_qpoint]
              #define NRML(i) row_n[i][i_qpoint]
              #define FLUX(i) flux(i_qpoint, i)
              FLUX(n_dim) = 0.;
              double pres = 0;
              for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                FLUX(n_dim) += READ(j_dim)*NRML(j_dim);
                pres += READ(j_dim)*READ(j_dim)/READ(n_dim);
              }
              pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
              double scaled = FLUX(n_dim)/READ(n_dim);
              FLUX(n_var - 1) = (READ(n_var - 1) + pres)*scaled;
              for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                FLUX(j_dim) = READ(j_dim)*scaled + pres*NRML(j_dim);
              }
              #undef FLUX
              #undef NRML
              #undef READ
            }

            // fetch boundary flux
            Eigen::Matrix<double, 2, n_var> boundary_values;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int is_positive : {0, 1}) {
                boundary_values(is_positive, i_var) = face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_qpoint/row_size + i_face_qpoint];
              }
            }

            Eigen::Matrix<double, row_size, n_var> row_w = -derivative(flux, boundary_values);

            // fetch jacobian determinant
            Eigen::Array<double, row_size, 1> row_det;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              row_det(i_qpoint) = determinant[i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }

            // Add dimensional component to update
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                int offset = i_outer*stride*row_size + i_inner + i_qpoint*stride;
                time_rate[i_var][offset] += row_w(i_qpoint, i_var)/row_det(i_qpoint);
              }
            }
            ++i_face_qpoint;
          }
        }
      }

      // write the updated solution
      for (int i_var = 0; i_var < n_var; ++i_var) {
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          const int i_dof = i_var*n_qpoint + i_qpoint;
          double updated = time_rate[i_var][i_qpoint]*d_t_by_d_pos*tss[i_qpoint] + state[i_dof];
          state[i_dof] = rkw*updated + (1. - rkw)*rk_reference[i_dof];
        }
      }
    }
  }
};

template<>
class Kernel_traits<Local_deformed>
{
  public:
  using base_t = Kernel<Deformed_element&>;
};

}
#endif
