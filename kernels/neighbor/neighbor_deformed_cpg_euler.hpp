#ifndef CARTDG_NEIGHBOR_DEFORMED_CPG_EULER_HPP_
#define CARTDG_NEIGHBOR_DEFORMED_CPG_EULER_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include "read_copy.hpp"
#define CARTDG_ATOMIC // avoids race condition since i_axis is not the same for all threads
#include "write_copy.hpp"
#undef CARTDG_ATOMIC
#include "hll_deformed_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void neighbor_deformed_cpg_euler(double** def_connections_r, double** def_connections_w,
                                 double** def_connections_j, int* i_axis, int* is_positive_face,
                                 int def_n_connections,
                                 Basis& basis, Kernel_settings& settings)
// "def_" (for "deformed") prepended to arguments to avoid naming conflicts in benchmark code
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  const int n_dim = n_var - 2;
  const int jac_size = n_dim*n_dim*n_face_qpoint;
  double mult = settings.d_t_by_d_pos/basis.node_weights()[0];
  double heat_rat = settings.cpg_heat_rat;

  #pragma omp parallel for
  for (int i_con = 0; i_con < def_n_connections; ++i_con)
  {
    double face_r [2*face_size];
    double face_jacobian [2*jac_size];
    double face_w [2*face_size];

    double** connect = def_connections_r + 2*i_con;
    for (int i_side : {0, 1})
    {
      int i_axis_side = i_axis[2*i_con + i_side];
      int stride = n_face_qpoint;
      for (int i = 0; i < i_axis_side; ++i) stride /= row_size;
      bool is_positive = is_positive_face[2*i_con + i_side] == 1;
      read_copy<n_var, n_qpoint, row_size>(connect[i_side], face_r + i_side*face_size, stride, is_positive);
      read_copy<n_dim*n_dim, n_qpoint, row_size>(def_connections_j[2*i_con + i_side], face_jacobian + i_side*jac_size, stride, is_positive);
    }

    if ((is_positive_face[2*i_con] != is_positive_face[2*i_con + 1]) && (i_axis[2*i_con] != i_axis[2*i_con + 1]))
    {
      if (n_dim == 3)
      {} // FIXME
      else
      {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var>> face_state (face_r);
        face_state.colwise().reverseInPlace();
        Eigen::Map<Eigen::Matrix<double, row_size, n_dim*n_dim>> face_jac (face_jacobian);
        face_jac.colwise().reverseInPlace();
      }
    }

    bool flip [] {is_positive_face[2*i_con] == 0, is_positive_face[2*i_con + 1] == 1};
    hll_deformed_cpg_euler<n_var - 2, n_face_qpoint>(face_r, face_w, face_jacobian, mult, i_axis + 2*i_con, flip, heat_rat);

    if ((is_positive_face[2*i_con] != is_positive_face[2*i_con + 1]) && (i_axis[2*i_con] != i_axis[2*i_con + 1]))
    {
      if (n_dim == 3)
      {}
      else
      {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var>> face_state (face_w);
        face_state.colwise().reverseInPlace();
      }
    }

    connect = def_connections_w + 2*i_con;

    for (int i_side : {0, 1})
    {
      int i_axis_side = i_axis[2*i_con + i_side];
      int stride = n_face_qpoint;
      for (int i = 0; i < i_axis_side; ++i) stride /= row_size;
      bool is_positive = is_positive_face[2*i_con + i_side];
      write_copy<n_var, n_qpoint, row_size>(face_w + i_side*face_size, connect[i_side], stride, is_positive);
    }
  }
}

}

#endif
