#ifndef CARTDG_NEIGHBOR_CONVECTIVE_HPP_
#define CARTDG_NEIGHBOR_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(regular, 3)
template<int n_var, int n_qpoint, int row_size>
void neighbor_convective(elem_con_vec& connections, Basis& basis, Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  double heat_rat = settings.cpg_heat_rat;

  for (unsigned stride = n_face_qpoint, i_dim = 0; stride > 0; stride /= row_size, ++i_dim)
  {
    #pragma omp parallel for
    for (unsigned i_con = 0; i_con < connections[i_dim].size(); ++i_con)
    {
      double face_r [2*face_size];
      double face_w [2*face_size];
      elem_con con = connections[i_dim][i_con];

      for (int i_side : {0, 1})
      {
        double* read = con[i_side]->face() + (i_dim*2 + 1 - i_side)*face_size;
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
        {
          face_r[i_side*face_size + i_face_dof] = read[i_face_dof];
        }
      }

      hll_cpg_euler<n_var - 2, n_face_qpoint>(face_r, face_w, 1., i_dim, heat_rat);

      for (int i_side : {0, 1})
      {
        double* write = con[i_side]->face() + (i_dim*2 + 1 - i_side)*face_size;
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
        {
          write[i_face_dof] = face_w[i_side*face_size + i_face_dof];
        }
      }
    }
  }
}

}
#endif
