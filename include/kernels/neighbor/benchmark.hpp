#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include "kernels/neighbor/read_copy.hpp"
#include "kernels/neighbor/write_copy.hpp"
#include "kernels/neighbor/average_flux.hpp"

template<n_var, int n_qpoint, int row_size>
void average_neighbor(double** connections_r, double** connections_w, int n_connections)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  for (int stride=n_face_qpoint, i_axis=0; stride > 0; stride /= row_size, ++i_axis)
  {
    double face0 [face_size];
    double face1 [face_size];
    double face_w [face_size];
    for (int i_con = 0; i_con < n_connections; ++i_con)
    {
      con_offset = i_axis*n_connections + 2*i_con;
      read_copy<n_var, n_qpoint, row_size>(connections_r[con_offset    ], face0, stride, 1);
      read_copy<n_var, n_qpoint, row_size>(connections_r[con_offset + 1], face1, stride, 0);
      average_flux<face_size>(&face0[0], &face1[0], &face_w[0]);
      write_copy<n_var, n_qpoint, row_size>(face_w, connections_w[con_offset    ], stride, 1);
      write_copy<n_var, n_qpoint, row_size>(face_w, connections_w[con_offset + 1], stride, 0);
    }
  }
}

#endif
