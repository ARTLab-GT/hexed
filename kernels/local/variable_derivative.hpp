#ifndef CARTDG_VARIABLE_DERIVATIVE_HPP_
#define CARTDG_VARIABLE_DERIVATIVE_HPP_

#include <vector>

#include <Eigen/Dense>

#include <math.hpp>
#include <Kernel_settings.hpp>
#include <Basis.hpp>

namespace cartdg
{

template<int n_dim, int row_size, bool add=false>
void variable_derivative(double* read, double* write, int i_dim,
                         Eigen::Matrix<double, row_size, row_size>& diff_mat, double d_pos)
{
  const int stride = custom_math::pow(row_size, n_dim - i_dim - 1);
  const int n_rows = custom_math::pow(row_size, i_dim);
  for (int i_outer = 0; i_outer < n_rows; ++i_outer)
  {
    for (int i_inner = 0; i_inner < stride; ++i_inner)
    {
      Eigen::Matrix<double, row_size, 1> row_r;
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
      {
        row_r(i_qpoint) = read[i_outer*stride*row_size + i_inner + i_qpoint*stride];
      }
      Eigen::Matrix<double, row_size, 1> row_w = diff_mat*row_r/d_pos;
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
      {
        int w_ind = i_outer*stride*row_size + i_inner + i_qpoint*stride;
        add ? write[w_ind] += row_w(i_qpoint) : write[w_ind] = row_w(i_qpoint);
      }
    }
  }
}

}

#endif
