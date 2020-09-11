#ifndef CPG_EULER_MATRIX_HPP_
#define CPG_EULER_MATRIX_HPP_

#include <Eigen/Dense>

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_matrix(double * read, double * write, int n_elem,
                      double * diff_mat, double * quad_weights,
                      double cfl, double sp_heat_rat = 1.4)
{

  // Fetch differentiation matrix
  double mat [row_size*row_size];
  for (int i_coef = 0; i_coef < row_size*row_size; ++i_coef)
  {
    mat[i_coef] = diff_mat[i_coef];
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

          // Calculate flux
          double flux [n_var*row_size];
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            #define READ(i) row_r[(i)*row_size + i_qpoint]
            #define FLUX(i) flux[(i)*row_size + i_qpoint]
            double veloc = READ(i_axis)/READ(n_var - 2);
            double pres = 0;
            for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
            {
              FLUX(j_axis) = READ(j_axis)*veloc;
              pres += READ(j_axis)*READ(j_axis)/READ(n_var - 2);
            }
            pres = (sp_heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
            FLUX(i_axis) += pres;
            FLUX(n_var - 2) = READ(i_axis);
            FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
            #undef FLUX
            #undef READ
          }

          // Differentiate flux
          Eigen::Map<Eigen::Matrix<double, row_size, n_var   >> f (&(flux[0]));
          Eigen::Map<Eigen::Matrix<double, row_size, n_var   >> w (&(row_w[0]));
          Eigen::Map<Eigen::Matrix<double, row_size, row_size>> d (&(mat[0]));
          w.noalias() = d*f*cfl;

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

#endif
