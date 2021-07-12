#ifndef CARTDG_LOCAL_DEFORMED_CPG_EULER_HPP_
#define CARTDG_LOCAL_DEFORMED_CPG_EULER_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_deformed_cpg_euler(double* read, double* write, double* jacobian, int n_elem, Basis& basis, Kernel_settings& settings)
{
  Eigen::Matrix<double, row_size, row_size> diff_mat = basis.diff_mat();
  double d_t_by_d_pos = settings.d_t_by_d_pos;
  double heat_rat = settings.cpg_heat_rat;
  const int n_dim = n_var - 2;

  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {

    double flux [n_dim][n_var][n_qpoint];
    double jacobian_determinant [n_qpoint];
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      // Get jacobian
      Eigen::Matrix<double, n_dim, n_dim> local_jacobian;
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
        {
          local_jacobian(i_axis, j_axis) = jacobian[((i_elem*n_dim + i_axis)*n_dim + j_axis)*n_qpoint + i_qpoint];
        }
      }
      jacobian_determinant[i_qpoint] = local_jacobian.determinant();

      // Compute flux in physical space
      double pres = 0;
      #define READ(i_var) read[(i_elem*n_var + i_var)*n_qpoint + i_qpoint]
      for (int i_axis = 0; i_axis < n_var - 2; ++i_axis)
      {
        pres += READ(i_axis)*READ(i_axis)/READ(n_var - 2);
      }
      pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);

      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        #define FLUX(i_var) flux[i_axis][i_var][i_qpoint]
        double veloc = READ(i_axis)/READ(n_var - 2);
        for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
        {
          FLUX(j_axis) = READ(j_axis)*veloc;
        }
        FLUX(i_axis) += pres;
        FLUX(n_var - 2) = READ(i_axis);
        FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
        #undef FLUX
      }
      #undef READ

      #if 0
      // Compute flux in reference space
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        double temp_flux [n_dim];
        for (int i_axis = 0; i_axis < n_dim; ++i_axis)
        {
          Eigen::Matrix<double, n_dim, n_dim> flux_mat = local_jacobian;
          for (int j_axis = 0; j_axis < n_dim; ++j_axis)
          {
            flux_mat(j_axis, i_axis) = flux[j_axis][i_var][i_qpoint];
          }
          temp_flux[i_axis] = flux_mat.determinant();
        }
        for (int i_axis = 0; i_axis < n_dim; ++i_axis)
        {
          flux[i_axis][i_var][i_qpoint] = temp_flux[i_axis];
        }
      }
      #endif
    }

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
          double phys_flux [n_dim][n_var][row_size];
          double row_flux [n_var][row_size];
          double row_jac [n_dim][n_dim][row_size];
          double row_w [n_var][row_size];
          for (int j_axis = 0; j_axis < n_dim; ++j_axis)
          {
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
              {
                phys_flux[j_axis][i_var][i_qpoint]
                = flux[j_axis][i_var][i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
          }
          for (int j_axis = 0; j_axis < n_dim; ++j_axis)
          {
            for (int k_axis = 0; k_axis < n_dim; ++k_axis)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
              {
                row_jac[j_axis][k_axis][i_qpoint] = jacobian[((i_elem*n_dim + j_axis)*n_dim + k_axis)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
          }

          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              Eigen::Matrix<double, n_dim, n_dim> flux_mat;
              for (int j_axis = 0; j_axis < n_dim; ++j_axis)
              {
                for (int k_axis = 0; k_axis < n_dim; ++k_axis)
                {
                  flux_mat(j_axis, k_axis) = row_jac[j_axis][k_axis][i_qpoint];
                }
              }
              for (int j_axis = 0; j_axis < n_dim; ++j_axis)
              {
                flux_mat(j_axis, i_axis) = phys_flux[j_axis][i_var][i_qpoint];
              }
              row_flux[i_var][i_qpoint] = flux_mat.determinant();
            }
          }

          // Differentiate flux
          Eigen::Map<Eigen::Matrix<double, row_size, n_var>> f (&(row_flux[0][0]));
          Eigen::Map<Eigen::Matrix<double, row_size, n_var>> w (&(row_w[0][0]));
          w.noalias() = -diff_mat*f*d_t_by_d_pos;

          // Write updated solution
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              const int i = + i_outer*stride*row_size + i_inner + i_qpoint*stride;
               write[(i_elem*n_var + i_var)*n_qpoint + i]
               += row_w[i_var][i_qpoint]/jacobian_determinant[i];
            }
          }

        }
      }
    }

  }
}

}

#endif
