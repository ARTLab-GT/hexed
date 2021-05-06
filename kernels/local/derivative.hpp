#ifndef CARTDG_DERIVATIVE_HPP_
#define CARTDG_DERIVATIVE_HPP_

#include <vector>

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>

namespace cartdg
{

template<int n_var_read, int n_var_write, int n_qpoint, int row_size, bool add=false>
void derivative(double* read, double* write, int n_elem, int i_var_read, int i_var_write, int i_axis,
                Basis& basis, Kernel_settings& settings)
{
  Eigen::Matrix<double, row_size, row_size> diff_mat = basis.diff_mat();

  int stride = n_qpoint/row_size;
  int n_rows = 1;
  for (int j_axis = 0; j_axis < i_axis; ++j_axis)
  {
    n_rows *= row_size;
    stride /= row_size;
  }

  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_outer = 0; i_outer < n_rows; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        // Fetch this row of data
        Eigen::Matrix<double, row_size, 1> row_r;
        for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
        {
          row_r(i_qpoint) = read[(i_elem*n_var_read + i_var_read)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
        }
        Eigen::Matrix<double, row_size, 1> row_w = diff_mat*row_r;
        for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
        {
          int w_ind = (i_elem*n_var_write + i_var_write)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride;
          add ? write[w_ind] += row_w(i_qpoint) : write[w_ind] = row_w(i_qpoint);
        }
      }
    }
  }
}

}

#endif
