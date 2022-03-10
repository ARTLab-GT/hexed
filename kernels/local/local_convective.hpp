#ifndef CARTDG_LOCAL_CONVECTIVE_HPP_
#define CARTDG_LOCAL_CONVECTIVE_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>
#include <Derivative.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(cartesian, 3)
template<int n_var, int n_qpoint, int row_size>
void local_convective(elem_vec& elements, Basis& basis, Kernel_settings& settings)
{
  Derivative<row_size> derivative (basis);
  const double d_t_by_d_pos = settings.d_t/settings.d_pos;
  const double heat_rat = settings.cpg_heat_rat;
  const int n_face_dof = n_var*n_qpoint/row_size;
  const double rk_weight = settings.rk_weight;

  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    double* state = elements[i_elem]->stage(0);
    double* rk_reference = state + n_var*n_qpoint;
    double time_rate [n_var][n_qpoint] {};
    double* face = elements[i_elem]->face();
    double* tss = elements[i_elem]->time_step_scale();

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
        state[i_dof] = rk_weight*updated + (1. - rk_weight)*rk_reference[i_dof];
      }
    }
  }
}

}

#endif
