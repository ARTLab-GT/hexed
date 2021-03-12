#ifndef CARTDG_CPG_EULER_MATRIX_HPP_
#define CARTDG_CPG_EULER_MATRIX_HPP_

#include <Eigen/Dense>

#include "../Kernel_settings.hpp"
#include "n_extrema.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_matrix(double * read, double * write, int n_elem,
                      Eigen::MatrixXd diff_mat_arg,
                      Eigen::VectorXd weights_arg,
                      Kernel_settings& settings)
{
  Eigen::Matrix<double, row_size, row_size> diff_mat = diff_mat_arg;
  double weights [row_size];
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) weights[i_qpoint] = weights_arg(i_qpoint);
  double d_t_by_d_pos = settings.d_t_by_d_pos;
  double heat_rat = settings.cpg_heat_rat;

  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {

    // Initialize updated solution to be equal to current solution
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      const int i = i_elem*n_qpoint*n_var + i_dof;
      write[i] = read[i];
    }

    // Perform update
    for (int stride = n_qpoint/row_size, n_rows = 1, i_axis = 0;
         n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_axis)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {

          // Fetch this row of data
          double row_r [n_var][row_size];
          double row_w [n_var][row_size];
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_r[i_var][i_qpoint]
              = read[(i_elem*n_var + i_var)*n_qpoint
                     + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }

          // Calculate flux
          double flux [n_var][row_size];
          double pressure [row_size];
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            #define READ(i) row_r[i][i_qpoint]
            #define FLUX(i) flux[i][i_qpoint]
            double veloc = READ(i_axis)/READ(n_var - 2);
            double& pres = pressure[i_qpoint];
            pres = 0;
            for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
            {
              FLUX(j_axis) = READ(j_axis)*veloc;
              pres += READ(j_axis)*READ(j_axis)/READ(n_var - 2);
            }
            pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
            FLUX(i_axis) += pres;
            FLUX(n_var - 2) = READ(i_axis);
            FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
            #undef FLUX
            #undef READ
          }

          // Differentiate flux
          Eigen::Map<Eigen::Matrix<double, row_size, n_var>> f (&(flux[0][0]));
          Eigen::Map<Eigen::Matrix<double, row_size, n_var>> w (&(row_w[0][0]));
          w.noalias() = -diff_mat*f*d_t_by_d_pos;

          // Extremum-counting limiter
          int n_ex = row_size;
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            n_ex = std::min<int>(n_ex, n_extrema<row_size>(row_w[i_var]));
          }
          if (n_ex > 1)
          {
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
              {
                row_w[i_var][i_qpoint] = 0.;
              }
              row_w[i_var][0] = flux[i_var][0]/weights[0]*d_t_by_d_pos;
              row_w[i_var][row_size - 1] = -flux[i_var][row_size - 1]/weights[row_size - 1]*d_t_by_d_pos;
            }
            double sound_speed = std::sqrt(heat_rat*pressure[0]/row_r[n_var - 2][0]);
            double veloc = row_r[i_axis][0]/row_r[n_var - 2][0];
            for (int i_qpoint = 1; i_qpoint < row_size; ++i_qpoint)
            {
              double wave_speed0 = veloc - sound_speed;
              double sound_speed = std::sqrt(heat_rat*pressure[i_qpoint]
                                              /row_r[n_var - 2][i_qpoint]);
              double veloc = row_r[i_axis][i_qpoint]/row_r[n_var - 2][i_qpoint];
              double wave_speed1 = veloc + sound_speed;
              for (int i_var = 0; i_var < n_var; ++i_var)
              {
                double num_flux;
                if      (std::min(wave_speed0, wave_speed1) >= 0) num_flux = flux[i_var][i_qpoint - 1];
                else if (std::max(wave_speed1, wave_speed0) <= 0) num_flux = flux[i_var][i_qpoint];
                else
                {
                  num_flux = (wave_speed1*flux[i_var][i_qpoint - 1] - wave_speed0*flux[i_var][i_qpoint]
                              + wave_speed0*wave_speed1*(row_r[i_var][i_qpoint] - row_r[i_var][i_qpoint - 1]))
                            / (wave_speed1 - wave_speed0);
                }
                row_w[i_var][i_qpoint - 1] -= num_flux/weights[i_qpoint - 1]*d_t_by_d_pos;
                row_w[i_var][i_qpoint    ] += num_flux/weights[i_qpoint    ]*d_t_by_d_pos;
              }
            }
          }

          // Write updated solution
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
               write[(i_elem*n_var + i_var)*n_qpoint
                     + i_outer*stride*row_size + i_inner + i_qpoint*stride]
               += row_w[i_var][i_qpoint];
            }
          }
        }
      }
    }

  }
}

}
#endif
