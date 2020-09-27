#ifndef CPG_EULER_COPY_HPP_
#define CPG_EULER_COPY_HPP_

#include <kernels/neighbor/read_copy.hpp>
#include <kernels/neighbor/write_copy.hpp>
#include <kernels/neighbor/cpg_euler_hll.hpp>

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_copy(double*** connections_r, double*** connections_w, int* n_connections,
                    Eigen::VectorXd weights_1d, double d_t_by_d_pos, double sp_heat_rat=1.4)
{
  const int n_face_qpoint = n_qpoint/row_size;
  double weights [n_qpoint];
  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) weights[i_qpoint] = 1.;
  for (int stride = 1, n_rows = 1; n_rows < n_qpoint;
       stride /= row_size, n_rows *= row_size)
  {
    for (int i_outer = 0; i_outer < n_rows; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
        {
          weights[(i_outer*row_size + i_qpoint)*stride + i_inner]
          *= weights_1d(i_qpoint)*d_t_by_d_pos;
        }
      }
    }
  }

  const int face_size = n_face_qpoint*n_var;
  for (int stride=n_face_qpoint, i_axis=0; stride > 0; stride /= row_size, ++i_axis)
  {
    double face_r [2*face_size];
    double face_w [2*face_size];
    for (int i_con = 0; i_con < n_connections[i_axis]; ++i_con)
    {
      double** connect = connections_r[i_axis];
      read_copy<n_var, n_qpoint, row_size>(connect[0], face_r            , stride, 1);
      read_copy<n_var, n_qpoint, row_size>(connect[1], face_r + face_size, stride, 0);

      cpg_euler_hll<n_var - 2, n_face_qpoint>(face_r, face_w, weights, i_axis, sp_heat_rat);

      connect = connections_w[i_axis];
      write_copy<n_var, n_qpoint, row_size>(face_w, connect[0], stride, 1);
      write_copy<n_var, n_qpoint, row_size>(face_w, connect[1], stride, 0);
    }
  }
}

#endif
