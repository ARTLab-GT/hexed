#ifndef CARTDG_JUMP_HPP_
#define CARTDG_JUMP_HPP_

#include "read_copy.hpp"
#include "write_copy.hpp"
#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var_read, int n_var_write, int n_qpoint, int row_size>
void jump(double** connections_r, double** connections_w, int n_con,
          int i_var_read, int i_var_write, int i_axis,
          const Eigen::VectorXd weights_1d, Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  double weight = weights_1d(0);
  int stride=n_face_qpoint;
  for (int j_axis = 0; j_axis < i_axis; ++j_axis) stride /= row_size;

  for (int i_con = 0; i_con < n_con; ++i_con)
  {
    double face_r [2][n_face_qpoint];
    double face_w [2][n_face_qpoint];

    double** connect = connections_r + 2*i_con;
    read_copy<1, n_qpoint, row_size>(connect[0] + i_var_read*n_qpoint, face_r[0], stride, 1);
    read_copy<1, n_qpoint, row_size>(connect[1] + i_var_read*n_qpoint, face_r[1], stride, 0);

    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
    {
      double avg = (face_r[0][i_qpoint] + face_r[1][i_qpoint])/2.;
      face_w[0][i_qpoint] =  (avg - face_r[0][i_qpoint])/weight;
      face_w[1][i_qpoint] = -(avg - face_r[1][i_qpoint])/weight;
    }

    connect = connections_w + 2*i_con;
    write_copy<1, n_qpoint, row_size>(face_w[0], connect[0] + i_var_write*n_qpoint, stride, 1);
    write_copy<1, n_qpoint, row_size>(face_w[1], connect[1] + i_var_write*n_qpoint, stride, 0);
  }
}

template<int n_var, int n_qpoint, int row_size>
void jump_r(double** connections_r, double** connections_w, int n_con,
           int i_var, int i_axis,
           const Eigen::VectorXd weights_1d, Kernel_settings& settings)
{
  jump<n_var, 1, n_qpoint, row_size>(connections_r, connections_w, n_con,
                                     i_var, 0, i_axis, weights_1d, settings);
}

template<int n_var, int n_qpoint, int row_size>
void jump_w(double** connections_r, double** connections_w, int n_con,
           int i_var, int i_axis,
           const Eigen::VectorXd weights_1d, Kernel_settings& settings)
{
  jump<1, n_var, n_qpoint, row_size>(connections_r, connections_w, n_con,
                                     0, i_var, i_axis, weights_1d, settings);
}

}

#endif
