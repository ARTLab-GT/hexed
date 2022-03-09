#ifndef CARTDG_LOCAL_GRADIENT_HPP_
#define CARTDG_LOCAL_GRADIENT_HPP_

#include <Basis.hpp>
#include <Element.hpp>
#include <Kernel_settings.hpp>

namespace cartdg
{

/*
 * Compute the gradient of state variable `i_var` (in stage 0) and write it to stage 2.
 * The `n_dim` components of the gradient are written to the first `n_dim` variables of the write storage.
 */
// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_gradient(elem_vec& elements, int i_var, Basis& basis, Kernel_settings& settings)
{
  const int n_face_dof = n_var*n_qpoint/row_size;
  // fetch basis properties
  const Eigen::Matrix<double, row_size, row_size> diff_mat {basis.diff_mat()};
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  Eigen::MatrixXd sign {{1, 0}, {0, -1}};
  const Eigen::Matrix<double, row_size, 1> inv_weights {Eigen::Array<double, row_size, 1>::Constant(1.)/basis.node_weights().array()};
  const Eigen::Matrix<double, row_size, 2> lift {inv_weights.asDiagonal()*basis.boundary().transpose()*sign};
  // fetch kernel parameters
  const double d_pos = settings.d_pos;
  // compute
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    double* read  = elements[i_elem]->stage(0);
    double* write = read + 2*n_var*n_qpoint;
    double* face = elements[i_elem]->face();
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      double* face0 = face + i_dim*2*n_face_dof;
      double* face1 = face0 + n_face_dof;
      int i_face_qpoint = 0;
      for (int i_outer = 0; i_outer < n_rows; ++i_outer) {
        for (int i_inner = 0; i_inner < stride; ++i_inner) {
          // Fetch this row of data
          Eigen::Matrix<double, row_size, 1> row_r;
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            row_r(i_qpoint) = read[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
          }
          // get neighbor correction
          const int face_offset = i_dim*n_qpoint/row_size + i_face_qpoint;
          Eigen::Matrix<double, 2, 1> boundary_values {face0[face_offset], face1[face_offset]};
          // compute derivative
          Eigen::Matrix<double, row_size, 1> row_w = diff_mat*row_r - lift*(boundary_values - boundary*row_r);
          // write row of data
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
            write[i_dim*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride] = row_w(i_qpoint)/d_pos;
          }
          ++i_face_qpoint;
        }
      }
    }
  }
}

}
#endif
