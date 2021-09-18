#ifndef CARTDG_CONT_VISC_HPP_
#define CARTDG_CONT_VISC_HPP_

#include <Eigen/Dense>

#include <math.hpp>
#include <Kernel_settings.hpp>
#include <Element.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void cont_visc(elem_con_vec& connections, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  constexpr int n_point = custom_math::pow(2, n_dim);
  constexpr int n_face_point = n_point/2;
  for (int stride = n_face_point, i_dim = 0; stride > 0; stride /= 2, ++i_dim)
  {
    #pragma omp parallel for
    for (unsigned i_con = 0; i_con < connections.size(); ++i_con)
    {
      double* visc [2] {connections[i_dim][i_con][0]->viscosity(),
                        connections[i_dim][i_con][1]->viscosity()};
      double face_r [2*n_face_point];
      double face_w [2*n_face_point];

      read_copy<1, n_point, 2>(visc[0], face_r               , stride, 1);
      read_copy<1, n_point, 2>(visc[1], face_r + n_face_point, stride, 0);

      for (int i_point = 0; i_point < n_face_point; ++i_point)
      {
        double max = std::max<double>(face_r[i_point], face_r[i_point + n_face_point]);
        face_w[i_point] = max - face_r[i_point];
        face_w[i_point + n_face_point] = max - face_r[i_point + n_face_point];
      }

      write_copy<1, n_point, 2>(face_w               , visc[0], stride, 1);
      write_copy<1, n_point, 2>(face_w + n_face_point, visc[1], stride, 0);
    }
  }
}

}
#endif
