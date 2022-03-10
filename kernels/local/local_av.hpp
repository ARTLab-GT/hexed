#ifndef CARTDG_LOCAL_GRADIENT_HPP_
#define CARTDG_LOCAL_GRADIENT_HPP_

#include <Basis.hpp>
#include <Element.hpp>
#include <Kernel_settings.hpp>
#include <Derivative.hpp>
#include <math.hpp>

namespace cartdg
{

/*
 * Compute the update to variable `i_var` due to artificial viscosity. Requires the
 * first `n_dim` variables of `stage(2)` to contain the components of the gradient.
 * Result is written to the `i_var`th variable of `stage(0)`.
 */
// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_av(elem_vec& elements, int i_var, Basis& basis, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  const int n_face_dof = n_var*n_qpoint/row_size;
  Derivative<row_size> derivative (basis);
  // fetch basis properties
  Eigen::Matrix<double, row_size, 2> interp;
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    interp(i_qpoint, 0) = 1. - basis.node(i_qpoint);
    interp(i_qpoint, 1) =      basis.node(i_qpoint);
  }
  Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  // fetch kernel parameters
  const double d_t_by_d_pos = settings.d_t/settings.d_pos;
  // compute
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    double* write = elements[i_elem]->stage(0);
    double* read  = write + 2*n_var*n_qpoint;
    double* face = elements[i_elem]->face();
    Eigen::Map<Eigen::Matrix<double, custom_math::pow(2, n_dim), 1>> vert_visc (elements[i_elem]->viscosity());
    Eigen::VectorXd visc = custom_math::hypercube_matvec(interp, vert_visc);
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      double* face0 = face + i_dim*2*n_face_dof;
      double* face1 = face0 + n_face_dof;
      int i_face_qpoint = 0;
      for (int i_outer = 0; i_outer < n_rows; ++i_outer) {
        for (int i_inner = 0; i_inner < stride; ++i_inner) {
          // Fetch this row of data
          Eigen::Matrix<double, row_size, 1> flux;
          Eigen::Matrix<double, row_size, 1> row_v;
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            int i_linear = i_outer*stride*row_size + i_inner + i_qpoint*stride;
            row_v(i_qpoint) = visc(i_linear);
            flux(i_qpoint) = read[i_dim*n_qpoint + i_linear]*row_v(i_qpoint);
          }
          // get neighbor correction
          const int face_offset = i_var*n_qpoint/row_size + i_face_qpoint;
          Eigen::Matrix<double, 2, 1> boundary_values {face0[face_offset], face1[face_offset]};
          boundary_values = boundary_values.cwiseProduct(boundary*row_v);
          // compute derivative
          Eigen::Matrix<double, row_size, 1> row_w = derivative(flux, boundary_values);
          // write row of data
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            write[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride] += row_w(i_qpoint)*d_t_by_d_pos;
          }
          ++i_face_qpoint;
        }
      }
    }
  }
}

}
#endif
