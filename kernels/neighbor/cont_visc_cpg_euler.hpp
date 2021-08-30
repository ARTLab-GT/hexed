#ifndef CARTDG_CONT_VISC_CPG_EULER_HPP_
#define CARTDG_CONT_VISC_CPG_EULER_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "../static_math.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void cont_visc_cpg_euler(double*** connections, int* n_connections, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  constexpr int n_point = static_math::pow(2, n_dim);
  constexpr int n_face_point = n_point/2;
  for (int stride = n_face_point, i_dim = 0; stride > 0; stride /= 2, ++i_dim)
  {
    #pragma omp parallel for
    for (int i_con = 0; i_con < n_connections[i_dim]; ++i_con)
    {
      double face_r [2*n_face_point];
      double face_w [2*n_face_point];
      double** connect = connections[i_dim] + 2*i_con;

      read_copy<1, n_point, 2>(connect[0], face_r               , stride, 1);
      read_copy<1, n_point, 2>(connect[1], face_r + n_face_point, stride, 0);

      for (int i_point = 0; i_point < n_face_point; ++i_point)
      {
        double max = std::max<double>(face_r[i_point], face_r[i_point + n_face_point]);
        face_w[i_point] = max - face_r[i_point];
        face_w[i_point + n_face_point] = max - face_r[i_point + n_face_point];
      }

      write_copy<1, n_point, 2>(face_w               , connect[0], stride, 1);
      write_copy<1, n_point, 2>(face_w + n_face_point, connect[1], stride, 0);
    }
  }
}

}
#endif
