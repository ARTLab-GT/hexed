#ifndef CARTDG_LOCAL_GRADIENT_HPP_
#define CARTDG_LOCAL_GRADIENT_HPP_

#include <Basis.hpp>
#include <Element.hpp>
#include <Kernel_settings.hpp>
#include <math.hpp>

namespace cartdg
{

/*
 * Compute the update to variable `i_var` due to artificial viscosity. Requires the
 * first `n_dim` variables of `stage(i_write)` to contain the components of the gradient.
 * Result is written to the `i_var`th variable of `stage(i_read)` (not `i_write`!).
 */
// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_av(elem_vec& elements, int i_var, Basis& basis, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  // fetch basis properties
  const Eigen::Matrix<double, row_size, row_size> diff_mat {basis.diff_mat()};
  Eigen::Matrix<double, row_size, 2> interp;
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    interp(i_qpoint, 0) = 1. - basis.node(i_qpoint);
    interp(i_qpoint, 1) =      basis.node(i_qpoint);
  }
  // fetch kernel parameters
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;
  const double d_t_by_d_pos = settings.d_t_by_d_pos;
  // compute
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    // going from `i_write` to `i_read`!
    double* read  = elements[i_elem]->stage(i_write);
    double* write = elements[i_elem]->stage(i_read);
    Eigen::Map<Eigen::Matrix<double, custom_math::pow(2, n_dim), 1>> vert_visc (elements[i_elem]->viscosity());
    Eigen::VectorXd visc = custom_math::hypercube_matvec(interp, vert_visc);
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer) {
        for (int i_inner = 0; i_inner < stride; ++i_inner) {
          // Fetch this row of data
          Eigen::Matrix<double, row_size, 1> flux;
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            int i_linear = i_outer*stride*row_size + i_inner + i_qpoint*stride;
            flux(i_qpoint) = read[i_dim*n_qpoint + i_linear]*visc(i_linear);
          }
          // compute derivative
          Eigen::Matrix<double, row_size, 1> row_w = diff_mat*flux;
          // write row of data
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            write[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride] += row_w(i_qpoint)*d_t_by_d_pos;
          }
        }
      }
    }
  }
}

}
#endif
