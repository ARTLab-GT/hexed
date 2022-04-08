#ifndef CARTDG_LOCAL_CARTESIAN_HPP_
#define CARTDG_LOCAL_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"

namespace cartdg
{

/*
 * Computes the local update for one Runge-Kutta stage.
 * Here, "local" means that all of the necessary data is contained within each `Element` object.
 * This includes the interior term, the face flux correction, and the Runge-Kutta update.
 * The numerical flux must already have been written to the face storage (which is the neighbor kernel's job).
 */
template <int n_dim, int row_size>
class Local_cartesian : public Kernel<Element&>
{
  Derivative<row_size> derivative;
  double dt;
  double rkw;
  const double heat_rat;

  public:
  Local_cartesian(const Basis& basis, double d_time, double rk_weight, double heat_ratio=1.4) :
    derivative{basis},
    dt{d_time},
    rkw{rk_weight},
    heat_rat{heat_ratio}
  {}

  virtual void operator()(Sequence<Element&>& elements)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    constexpr int n_face_dof = n_var*n_qpoint/row_size;

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      double* state = elements[i_elem].stage(0);
      double* rk_reference = state + n_var*n_qpoint;
      double time_rate [n_var][n_qpoint] {};
      double* face = elements[i_elem].face();
      double* tss = elements[i_elem].time_step_scale();
      const double d_t_by_d_pos = dt/elements[i_elem].nominal_size();

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
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                row_r[i_var][i_qpoint] = state[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
            // fetch boundary flux
            Eigen::Matrix<double, 2, n_var> boundary_values;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int is_positive : {0, 1}) {
                boundary_values(is_positive, i_var) = face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_qpoint/row_size + i_face_qpoint];
              }
            }

            // Calculate flux
            Eigen::Matrix<double, row_size, n_var> flux;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              #define READ(i) row_r[i][i_qpoint]
              #define FLUX(i) flux(i_qpoint, i)
              double veloc = READ(i_dim)/READ(n_var - 2);
              double pres = 0;
              for (int j_dim = 0; j_dim < n_var - 2; ++j_dim) {
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
            // Differentiate flux
            Eigen::Matrix<double, row_size, n_var> row_w = -derivative(flux, boundary_values);

            // Add dimensional component to update
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                int offset = i_outer*stride*row_size + i_inner + i_qpoint*stride;
                time_rate[i_var][offset] += row_w(i_qpoint, i_var);
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
class Kernel_traits<Local_cartesian>
{
  public:
  using base_t = Kernel<Element&>;
};

}
#endif