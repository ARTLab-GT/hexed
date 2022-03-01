#ifndef CARTDG_LOCAL_GRADIENT_HPP_
#define CARTDG_LOCAL_GRADIENT_HPP_

#include <Basis.hpp>
#include <Element.hpp>
#include <Kernel_settings.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_gradient(elem_vec& elements, int i_var, Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, row_size, row_size> diff_mat {basis.diff_mat()};
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    double* read  = elements[i_elem]->stage(i_read);
    double* write = elements[i_elem]->stage(i_write);
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer) {
        for (int i_inner = 0; i_inner < stride; ++i_inner) {
          // Fetch this row of data
          Eigen::Matrix<double, row_size, 1> row_r;
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            row_r(i_qpoint) = read[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
          }
          // compute derivative
          Eigen::Matrix<double, row_size, 1> row_w = diff_mat*row_r;
          // write row of data
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            write[i_dim*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride] = row_w(i_qpoint);
          }
        }
      }
    }
  }
}

}
#endif
