#ifndef CARTDG_LOCAL_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_LOCAL_DEFORMED_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void local_deformed_convective(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  Eigen::Matrix<double, row_size, row_size> diff_mat = basis.diff_mat();
  double d_t_by_d_pos = settings.d_t_by_d_pos;
  double heat_rat = settings.cpg_heat_rat;
  const int n_dim = n_var - 2;
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;

  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < def_elements.size(); ++i_elem)
  {
    double* read     = def_elements[i_elem]->stage(i_read);
    double* write    = def_elements[i_elem]->stage(i_write);
    double* jacobian = def_elements[i_elem]->jacobian();

    // Precompute flux for all qpoints
    double flux [n_dim][n_var][n_qpoint];
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      // Compute flux in physical space
      double pres = 0;
      #define READ(i_var) read[(i_var)*n_qpoint + i_qpoint]
      for (int i_dim = 0; i_dim < n_var - 2; ++i_dim)
      {
        pres += READ(i_dim)*READ(i_dim)/READ(n_var - 2);
      }
      pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);

      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        #define FLUX(i_var) flux[i_dim][i_var][i_qpoint]
        double veloc = READ(i_dim)/READ(n_var - 2);
        for (int j_dim = 0; j_dim < n_var - 2; ++j_dim)
        {
          FLUX(j_dim) = READ(j_dim)*veloc;
        }
        FLUX(i_dim) += pres;
        FLUX(n_var - 2) = READ(i_dim);
        FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
        #undef FLUX
      }
      #undef READ
    }

    // Initialize updated solution to be equal to current solution
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      write[i_dof] = read[i_dof];
    }

    // Perform update
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0;
         n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {

          // Fetch this row of data
          double phys_flux [n_dim][n_var][row_size];
          double row_jac [n_dim][n_dim][row_size];
          double row_w [n_var][row_size];
          for (int j_dim = 0; j_dim < n_dim; ++j_dim)
          {
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
              {
                phys_flux[j_dim][i_var][i_qpoint]
                = flux[j_dim][i_var][i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
          }
          for (int j_dim = 0; j_dim < n_dim; ++j_dim)
          {
            for (int k_dim = 0; k_dim < n_dim; ++k_dim)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
              {
                row_jac[j_dim][k_dim][i_qpoint] = jacobian[(j_dim*n_dim + k_dim)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
          }

          // differentiate flux
          double flux_diff [n_dim][n_var][row_size];
          Eigen::Map<Eigen::Matrix<double, row_size, n_dim*n_var>> pf (&(phys_flux[0][0][0]));
          Eigen::Map<Eigen::Matrix<double, row_size, n_dim*n_var>> fd (&(flux_diff[0][0][0]));
          fd.noalias() = diff_mat*pf;

          double jac_det [row_size];
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            Eigen::Matrix<double, n_dim, n_dim> flux_mat;
            for (int j_dim = 0; j_dim < n_dim; ++j_dim)
            {
              for (int k_dim = 0; k_dim < n_dim; ++k_dim)
              {
                flux_mat(j_dim, k_dim) = row_jac[j_dim][k_dim][i_qpoint];
              }
            }
            jac_det[i_qpoint] = flux_mat.determinant();
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int j_dim = 0; j_dim < n_dim; ++j_dim)
              {
                flux_mat(j_dim, i_dim) = flux_diff[j_dim][i_var][i_qpoint];
              }
              row_w[i_var][i_qpoint] = flux_mat.determinant()*-d_t_by_d_pos;
            }
          }

          // Write updated solution
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              const int i = + i_outer*stride*row_size + i_inner + i_qpoint*stride;
              write[i_var*n_qpoint + i] += row_w[i_var][i_qpoint]/jac_det[i_qpoint];
            }
          }
        }
      }
    }
  }
}

}

#endif
