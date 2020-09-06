#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include <Eigen/Dense>

template<int n_var, int n_qpoint, int row_size>
void copy(double * read, double * write, int n_elem,
          double * diff_mat)
{
  for (int i_elem = 0; i_elem < n_elem*n_var; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      write[i_elem*n_qpoint + i_qpoint] = 0.;
    }
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      write[i_elem*n_qpoint + i_qpoint] += read[i_elem*n_qpoint + i_qpoint];
    }
  }
}

template<int n_var, int n_qpoint, int row_size>
void basic_tensor(double * read, double * write, int n_elem,
                  double * diff_mat)
{
  for (int i_elem = 0; i_elem < n_elem*n_var; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      write[i_elem*n_qpoint + i_qpoint] = 0.;
    }
    for (int stride = 1, n_rows = n_qpoint/row_size;
         stride < n_qpoint;
         stride *= row_size, n_rows /= row_size)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            for (int j_qpoint = 0; j_qpoint < row_size; ++j_qpoint)
            {
              write  [i_elem*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride]
              += read[i_elem*n_qpoint + i_outer*stride*row_size + i_inner + j_qpoint*stride];
            }
          }
        }
      }
    }
  }
}

template<int n_var, int n_qpoint, int row_size>
void update_add(double * read, double * write, int n_elem,
                double * diff_mat)
{
  double mat [row_size*row_size];
  for (int i_coef = 0; i_coef < row_size*row_size; ++i_coef)
  {
    mat[i_coef] = diff_mat[i_coef];
  }
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      const int i = i_elem*n_qpoint*n_var + i_dof;
      write[i] = read[i];
    }
    for (int stride = 1, n_rows = n_qpoint/row_size, i_axis = 0;
         stride < n_qpoint;
         stride *= row_size, n_rows /= row_size, ++i_axis)
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
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_w[i_var*row_size + i_qpoint] = 0.;
              for (int j_qpoint = 0; j_qpoint < row_size; ++j_qpoint)
              {
                row_w[i_var*row_size + i_qpoint] += row_r[i_var*row_size + j_qpoint];
              }
            }
          }
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

template<int n_var, int n_qpoint, int row_size>
void update_matvec(double * read, double * write, int n_elem,
                   double * diff_mat)
{
  double mat [row_size*row_size];
  for (int i_coef = 0; i_coef < row_size*row_size; ++i_coef)
  {
    mat[i_coef] = diff_mat[i_coef];
  }
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      const int i = i_elem*n_qpoint*n_var + i_dof;
      write[i] = read[i];
    }
    for (int stride = 1, n_rows = n_qpoint/row_size, i_axis = 0;
         stride < n_qpoint;
         stride *= row_size, n_rows /= row_size, ++i_axis)
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
          Eigen::Map<Eigen::Matrix<double, row_size, 1>> r (&(row_r[0]));
          Eigen::Map<Eigen::Matrix<double, row_size, 1>> w (&(row_w[0]));
          Eigen::Map<Eigen::Matrix<double, row_size, row_size>> dm (&(mat[0]));
          w.noalias() = dm*r;
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
