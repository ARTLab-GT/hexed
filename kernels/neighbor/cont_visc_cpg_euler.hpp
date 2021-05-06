#ifndef CARTDG_CONT_VISC_CPG_EULER_HPP_
#define CARTDG_CONT_VISC_CPG_EULER_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "../static_math.hpp"

namespace cartdg
{

// AUTOGENERATE
template<int n_var, int n_qpoint, int row_size>
void cont_visc_cpg_euler(double*** connections, int* n_connections, Kernel_settings& settings)
{
  #if 0
  const int n_face_qpoint = n_qpoint/row_size;
  double mult = settings.d_t_by_d_pos/weights_1d(0);
  double heat_rat = settings.cpg_heat_rat;

  const int face_size = n_face_qpoint*n_var;
  for (int stride=n_face_qpoint, i_axis=0; stride > 0; stride /= row_size, ++i_axis)
  {
    #pragma omp parallel for
    for (int i_con = 0; i_con < n_connections[i_axis]; ++i_con)
    {
      double face_r [2*face_size];
      double face_w [2*face_size];

      double** connect = connections_r[i_axis] + 2*i_con;
      read_copy<n_var, n_qpoint, row_size>(connect[0], face_r            , stride, 1);
      read_copy<n_var, n_qpoint, row_size>(connect[1], face_r + face_size, stride, 0);

      cpg_euler_hll<n_var - 2, n_face_qpoint>(face_r, face_w, mult, i_axis, heat_rat);

      connect = connections_w[i_axis] + 2*i_con;
      write_copy<n_var, n_qpoint, row_size>(face_w            , connect[0], stride, 1);
      write_copy<n_var, n_qpoint, row_size>(face_w + face_size, connect[1], stride, 0);
    }
  }
  #endif
}

}
#endif
