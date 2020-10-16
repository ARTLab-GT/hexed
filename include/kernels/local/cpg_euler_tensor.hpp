#ifndef CARTDG_CPG_EULER_TENSOR_HPP_
#define CARTDG_CPG_EULER_TENSOR_HPP_

#include <Eigen/Dense>

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_tensor(double * read, double * write, int n_elem,
                      const Eigen::MatrixXd& diff_mat, double d_t_by_d_x, double sp_heat_rat = 1.4)
{

  // Fetch differentiation matrix
  double mat [row_size*row_size];
  for (int row = 0; row < row_size; ++row)
  {
    for (int col = 0; col < row_size; ++col)
    {
      mat[row*row_size + col] = diff_mat(row, col);
    }
  }

  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {

    // Initialize updated solution to be equal to current solution
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      const int i = i_elem*n_qpoint*n_var + i_dof;
      write[i] = read[i];
    }

    // Perform update
    for (int stride = 1, n_rows = n_qpoint/row_size, i_axis = 0;
         stride < n_qpoint;
         stride *= row_size, n_rows /= row_size, ++i_axis)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {

          // Fetch this row of data
          double row_r [row_size*n_var];
          double row_w [row_size*n_var];
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_r[row_size*i_var + i_qpoint]
              = read[(i_elem*n_var + i_var)*n_qpoint
                     + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }

          // Calculate intermediary variables related to the flux
          double flux_var [(n_var + 2)*row_size];
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            #define READ(i) row_r[(i)*row_size + i_qpoint]
            #define FLUX(i) flux_var[(i)*row_size + i_qpoint]
            FLUX(n_var - 2) = READ(i_axis);
            FLUX(n_var - 1) = 0;
            for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
            {
              FLUX(j_axis) = READ(j_axis)/READ(n_var - 2);
              FLUX(n_var - 1) += FLUX(j_axis)*READ(j_axis);
            }
            FLUX(n_var - 1) = READ(i_axis)/2*(sp_heat_rat - 1.)*(READ(n_var - 1) - 0.5*FLUX(n_var - 1));
            FLUX(n_var) = std::log(READ(n_var - 2));
            FLUX(n_var + 1) = std::log(READ(n_var - 1));
            #undef FLUX
            #undef READ
          }

          // Differentiate flux
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
             row_w[i_var*row_size + i_qpoint] = 0;
            }
            for (int j_qpoint = 0; j_qpoint < row_size; ++j_qpoint)
            {
              #define WRITE(i) row_w[(i)*row_size + i_qpoint]
              double coef = d_t_by_d_x*mat[i_qpoint*row_size + j_qpoint];
              double mass_log = (flux_var[(n_var - 2)*row_size + i_qpoint] - flux_var[(n_var - 2)*row_size + j_qpoint])/
                                (flux_var[(n_var - 0)*row_size + i_qpoint] - flux_var[(n_var - 0)*row_size + j_qpoint]);
              double beta_log = (flux_var[(n_var - 1)*row_size + i_qpoint] - flux_var[(n_var - 1)*row_size + j_qpoint])/
                                (flux_var[(n_var + 1)*row_size + i_qpoint] - flux_var[(n_var + 1)*row_size + j_qpoint]);
              double mass_ari = 0.5*(flux_var[(n_var - 2)*row_size + i_qpoint] + flux_var[(n_var - 2)*row_size + j_qpoint]);
              double beta_ari = 0.5*(flux_var[(n_var - 1)*row_size + i_qpoint] + flux_var[(n_var - 1)*row_size + j_qpoint]);
              double veloc_ari = 0.5*(flux_var[i_axis*row_size + i_qpoint] + flux_var[i_axis*row_size + j_qpoint]);
              WRITE(n_var - 2) += coef*mass_log*veloc_ari;
              double kin_ener = 0;
              for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
              {
                double veloc_ari_j = 0.5*(flux_var[j_axis*row_size + i_qpoint] + flux_var[j_axis*row_size + j_qpoint]);
                WRITE(j_axis) = WRITE(n_var - 2)*veloc_ari_j + mass_ari/(2*beta_ari);
                kin_ener += 0.5*(flux_var[j_axis*row_size + i_qpoint]*flux_var[j_axis*row_size + i_qpoint] + flux_var[j_axis*row_size + j_qpoint]*flux_var[j_axis*row_size + j_qpoint]);
              }
              WRITE(n_var - 1) = (0.5/((sp_heat_rat - 1.)*beta_log) - 0.5*kin_ener)*WRITE(n_var - 2) + 0;
            }
          }

          // Write updated solution
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
               write[(i_elem*n_var + i_var)*n_qpoint
                     + i_outer*stride*row_size + i_inner + i_qpoint*stride]
               += row_w[row_size*i_var + i_qpoint];
            }
          }

        }
      }
    }

  }
}

}
#endif
