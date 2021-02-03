#ifndef CARTDG_CPG_EULER_COPY_HPP_
#define CARTDG_CPG_EULER_COPY_HPP_

#include <kernels/neighbor/read_copy.hpp>
#include <kernels/neighbor/write_copy.hpp>
#include <kernels/neighbor/cpg_euler_hll.hpp>

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_copy(double*** connections_r, double*** connections_w, int* n_connections,
                    const Eigen::VectorXd weights_1d,
                    double d_t_by_d_pos, double sp_heat_rat=1.4)
{
  const int n_face_qpoint = n_qpoint/row_size;
  double mult = d_t_by_d_pos/weights_1d(0);

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

      cpg_euler_hll<n_var - 2, n_face_qpoint>(face_r, face_w, mult, i_axis, sp_heat_rat);

      connect = connections_w[i_axis] + 2*i_con;
      write_copy<n_var, n_qpoint, row_size>(face_w            , connect[0], stride, 1);
      write_copy<n_var, n_qpoint, row_size>(face_w + face_size, connect[1], stride, 0);
    }
  }
}

}
#endif
