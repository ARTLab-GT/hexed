#ifndef CARTDG_CPG_EULER_COPY_DEFORMED_HPP_
#define CARTDG_CPG_EULER_COPY_DEFORMED_HPP_

#include "read_copy.hpp"
#include "write_copy.hpp"
#include "cpg_euler_hll_deformed.hpp"
#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_copy_deformed(double** connections_r, double** connections_w,
                             double** jacobian, int* i_axis, int* is_positive_face,
                             int n_connections,
                             const Eigen::VectorXd weights_1d, Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  const int n_dim = n_var - 2;
  const int jac_size = n_dim*n_dim*n_face_qpoint;
  double mult = settings.d_t_by_d_pos/weights_1d(0);
  double heat_rat = settings.cpg_heat_rat;

  // FIXME: fix race condition to allow parallelism
  //#pragma omp parallel for
  for (int i_con = 0; i_con < n_connections; ++i_con)
  {
    double face_r [2*face_size];
    double face_jacobian[2*jac_size];
    double face_w [2*face_size];

    double** connect = connections_r + 2*i_con;
    for (int i_side : {0, 1})
    {
      int i_axis_side = i_axis[2*i_con + i_side];
      int stride = n_face_qpoint;
      for (int i = 0; i < i_axis_side; ++i) stride /= row_size;
      bool is_positive = is_positive_face[2*i_con + i_side] == 1;
      read_copy<n_var, n_qpoint, row_size>(connect[i_side], face_r + i_side*face_size, stride, is_positive);
      read_copy<n_dim*n_dim, n_qpoint, row_size>(jacobian[2*i_con + i_side], face_jacobian + i_side*jac_size, stride, is_positive);
    }

    if ((is_positive_face[2*i_con] != is_positive_face[2*i_con + 1]) && (i_axis[2*i_con] != i_axis[2*i_con + 1]))
    {
      if (n_dim == 3)
      {}
      else
      {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var>> face_state (face_r);
        face_state.colwise().reverseInPlace();
        Eigen::Map<Eigen::Matrix<double, row_size, n_dim*n_dim>> face_jac (face_jacobian);
        face_jac.colwise().reverseInPlace();
      }
    }

    bool flip [] {is_positive_face[2*i_con] == 0, is_positive_face[2*i_con + 1] == 1};
    cpg_euler_hll_deformed<n_var - 2, n_face_qpoint>(face_r, face_w, face_jacobian, mult, i_axis + 2*i_con, flip, heat_rat);

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

    connect = connections_w + 2*i_con;

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