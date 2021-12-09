#ifndef CARTDG_BENCHMARK_HPP_
#define CARTDG_BENCHMARK_HPP_

#include <Eigen/Dense>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK
template<int n_var, int n_qpoint, int row_size>
void copy(double* read, double* write, int n_elem)
{
  const int n_dof = n_var*n_qpoint;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    for (int i_dof = 0; i_dof < n_dof; ++i_dof) {
      write[i_elem*n_dof + i_dof] = read[i_elem*n_dof + i_dof];
    }
  }
}

template<int n_var, int n_qpoint, int row_size>
void update_matvec(double * read, double * write, int n_elem,
                   const Eigen::MatrixXd& diff_mat_arg)
{
  Eigen::Matrix<double, row_size, row_size> diff_mat = diff_mat_arg;

  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      const int i = i_elem*n_qpoint*n_var + i_dof;
      write[i] = read[i];
    }
    for (int stride = 1, n_rows = n_qpoint/row_size, i_dim = 0;
         stride < n_qpoint;
         stride *= row_size, n_rows /= row_size, ++i_dim)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
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
          Eigen::Map<Eigen::Matrix<double, row_size, n_var   >> r  (&(row_r[0]));
          Eigen::Map<Eigen::Matrix<double, row_size, n_var   >> w  (&(row_w[0]));
          w.noalias() = diff_mat*r;
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
