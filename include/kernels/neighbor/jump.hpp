#ifndef CARTDG_JUMP_HPP_
#define CARTDG_JUMP_HPP_

#include "read_copy.hpp"
#include "write_copy.hpp"
#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void jump(std::vector<int>& con_inds, double** connections_r, double** connections_w,
          int i_var, int i_axis,
          const Eigen::VectorXd weights_1d, Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  double weight = weights_1d(0);
  int stride=n_face_qpoint;
  for (int j_axis = 0; j_axis < i_axis; ++j_axis) stride /= row_size;

  for (int i_con : con_inds)
  {
    double face_r [2][n_face_qpoint];
    double face_w [2][n_face_qpoint];

    double** connect = connections_r + 2*i_con;
    read_copy<1, n_qpoint, row_size>(connect[0] + i_var*n_qpoint, face_r[0], stride, 1);
    read_copy<1, n_qpoint, row_size>(connect[1] + i_var*n_qpoint, face_r[1], stride, 0);

    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
    {
      double avg = (face_r[0][i_qpoint] + face_r[1][i_qpoint])/2.;
      face_w[0][i_qpoint] =  (avg - face_r[0][i_qpoint])/weight;
      face_w[1][i_qpoint] = -(avg - face_r[1][i_qpoint])/weight;
    }

    connect = connections_w + 2*i_con;
    write_copy<1, n_qpoint, row_size>(face_w[0], connect[0], stride, 1);
    write_copy<1, n_qpoint, row_size>(face_w[1], connect[1], stride, 0);
  }
}

}

#endif
