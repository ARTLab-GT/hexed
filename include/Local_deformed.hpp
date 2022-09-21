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
      double flux [n_dim][n_var][n_qpoint];
      const double d_t_by_d_pos = dt/elem.nominal_size();

      // precompute flux
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          #define READ(i) state[(i)*n_qpoint + i_qpoint]
          #define FLUX(i) flux[i_dim][i][i_qpoint]
          double veloc = READ(i_dim)/READ(n_var - 2);
          double pres = 0;
          for (int j_dim = 0; j_dim < n_var - 2; ++j_dim)
          {
            FLUX(j_dim) = READ(j_dim)*veloc;
            pres += READ(j_dim)*READ(j_dim)/READ(n_var - 2);
          }
          pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
          FLUX(i_dim) += pres;
          FLUX(n_var - 2) = READ(i_dim);
          FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
          #undef FLUX
          #undef READ
        }
      }

      // Compute update
      for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
           stride /= row_size, n_rows *= row_size, ++i_dim)
      {
        int i_face_qpoint {0};
        for (int i_outer = 0; i_outer < n_rows; ++i_outer)
        {
          for (int i_inner = 0; i_inner < stride; ++i_inner)
          {
            // fetch flux
            Eigen::Matrix<double, row_size, n_var> row_f;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              for (int i_var = 0; i_var < n_var; ++i_var) {
                double f = 0.;
                for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                  f += flux[j_dim][i_var][i_outer*stride*row_size + i_inner + i_qpoint*stride]
                       *normals[(i_dim*n_dim + j_dim)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
                }
                row_f(i_qpoint, i_var) = f;
              }
            }

            // fetch jacobian determinant
            Eigen::Array<double, row_size, 1> row_det;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              row_det(i_qpoint) = determinant[i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }

            // fetch boundary flux
            Eigen::Matrix<double, 2, n_var> boundary_values;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int is_positive : {0, 1}) {
                boundary_values(is_positive, i_var) = face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_qpoint/row_size + i_face_qpoint];
              }
            }

            Eigen::Matrix<double, row_size, n_var> row_w = -derivative(row_f, boundary_values);

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
