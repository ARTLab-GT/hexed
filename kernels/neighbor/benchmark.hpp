#ifndef CARTDG_BENCHMARK_HPP_
#define CARTDG_BENCHMARK_HPP_

#include <Eigen/Dense>

#include <Basis.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "average.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void average_neighbor(double*** connections_r, double*** connections_w, int* n_connections,
                      Basis& basis, double d_t_by_d_pos)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  for (int stride = n_face_qpoint, i_dim = 0; stride > 0; stride /= row_size, ++i_dim)
  {
    double face0 [face_size];
    double face1 [face_size];
    double face_w [face_size];

    #pragma omp parallel for
    for (int i_con = 0; i_con < n_connections[i_dim]; ++i_con)
    {
      read_copy<n_var, n_qpoint, row_size>(connections_r[i_dim][2*i_con], face0, stride, 1);
      read_copy<n_var, n_qpoint, row_size>(connections_r[i_dim][2*i_con + 1], face1, stride, 0);
      average<face_size>(&face0[0], &face1[0], &face_w[0]);
      write_copy<n_var, n_qpoint, row_size>(face_w, connections_w[i_dim][2*i_con], stride, 1);
      write_copy<n_var, n_qpoint, row_size>(face_w, connections_w[i_dim][2*i_con + 1],
                                            stride, 0);
    }
  }
}

}
#endif
