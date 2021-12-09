#ifndef CARTDG_NEIGHBOR_CONVECTIVE_HPP_
#define CARTDG_NEIGHBOR_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(cartesian, 3)
template<int n_var, int n_qpoint, int row_size>
void neighbor_convective(elem_con_vec& connections, Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  double heat_rat = settings.cpg_heat_rat;

  for (unsigned stride = n_face_qpoint, i_dim = 0; stride > 0; stride /= row_size, ++i_dim)
  {
    #pragma omp parallel for
    for (unsigned i_con = 0; i_con < connections[i_dim].size(); ++i_con)
    {
      double face [2*face_size];
      elem_con con = connections[i_dim][i_con];

      for (int i_side : {0, 1})
      {
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
        {
          face[i_side*face_size + i_face_dof] = con[i_side][i_face_dof];
        }
      }

      hll_cpg_euler<n_var - 2, n_face_qpoint>(face, face, 1., i_dim, heat_rat);

      for (int i_side : {0, 1})
      {
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
        {
          con[i_side][i_face_dof] = face[i_side*face_size + i_face_dof];
        }
      }
    }
  }
}

}
#endif
