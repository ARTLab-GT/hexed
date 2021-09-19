#ifndef CARTDG_VARIABLE_JUMP_HPP_
#define CARTDG_VARIABLE_JUMP_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <math.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

template<int n_qpoint, int row_size>
void variable_jump(std::array<double*, 2> read, std::array<double*, 2> write,
                   int i_dim, double weight)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int stride = n_face_qpoint/custom_math::pow(row_size, i_dim);

  double face_r [2][n_face_qpoint] {}; // initialize to avoid "-Wmaybe-uninitialized"
  double face_w [2][n_face_qpoint] {};

  read_copy<1, n_qpoint, row_size>(read[0], face_r[0], stride, 1);
  read_copy<1, n_qpoint, row_size>(read[1], face_r[1], stride, 0);

  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double avg = (face_r[0][i_qpoint] + face_r[1][i_qpoint])/2.;
    face_w[0][i_qpoint] =  (avg - face_r[0][i_qpoint])/weight;
    face_w[1][i_qpoint] = -(avg - face_r[1][i_qpoint])/weight;
  }

  write_copy<1, n_qpoint, row_size>(face_w[0], write[0], stride, 1);
  write_copy<1, n_qpoint, row_size>(face_w[1], write[1], stride, 0);
}

}

#endif
