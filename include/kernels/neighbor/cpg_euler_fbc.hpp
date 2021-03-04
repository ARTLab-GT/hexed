#ifndef CARTDG_CPG_EULER_FBC_HPP_
#define CARTDG_CPG_EULER_FBC_HPP_

#include "../Kernel_settings.hpp"
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

template <int n_var, int n_qpoint, int row_size>
void cpg_euler_fbc(double* read, double* write, double* jacobian, int* i_elem,
                   int* i_dim, int* is_positive_face, int n_bc, double weight,
                   Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  const int n_dim = n_var - 2;
  const int jac_size = n_dim*n_dim*n_face_qpoint;
  double mult = settings.d_t_by_d_pos/weight;
  double heat_rat = settings.cpg_heat_rat;

  for (int i_bc = 0; i_bc < n_bc; ++i_bc)
  {
    int stride = n_face_qpoint;
    for (int i = 0; i < i_dim[i_bc]; ++i) stride /= row_size;

    double face_r [face_size];
    double face_jacobian[jac_size];
    double face_w [face_size] {};
    read_copy<n_var, n_qpoint, row_size>(read + n_var*n_qpoint*i_bc, face_r, stride,
                                         is_positive_face[i_bc]);
    read_copy<n_dim*n_dim, n_qpoint, row_size>(jacobian + n_dim*n_dim*n_qpoint*i_bc, face_jacobian, stride, is_positive_face[i_bc]);

    write_copy<n_var, n_qpoint, row_size>(face_w, write + n_var*n_qpoint*i_bc, stride,
                                          is_positive_face[i_bc]);
  }
}

}

#endif