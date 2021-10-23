#ifndef CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void neighbor_deformed_convective(def_elem_con_vec& def_connections, Basis& basis, Kernel_settings& settings)
// "def_" (for "deformed") prepended to arguments to avoid naming conflicts in benchmark code
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  double heat_rat = settings.cpg_heat_rat;

  #pragma omp parallel for
  for (unsigned i_con = 0; i_con < def_connections.size(); ++i_con)
  {
    double face_r [2*face_size];
    double face_w [2*face_size];
    Deformed_elem_con con = def_connections[i_con];

    for (int i_side : {0, 1})
    {
      double* read = con.element[i_side]->face() + (con.i_dim[i_side]*2 + 1 - i_side)*face_size;
      for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
      {
        face_r[i_side*face_size + i_face_dof] = read[i_face_dof];
      }
    }

    hll_cpg_euler<n_var - 2, n_face_qpoint>(face_r, face_w, 1., con.i_dim[0], heat_rat);

    for (int i_side : {0, 1})
    {
      double* write = con.element[i_side]->face() + (con.i_dim[i_side]*2 + 1 - i_side)*face_size;
      for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
      {
        write[i_face_dof] = face_w[i_side*face_size + i_face_dof];
      }
    }
  }
}

}

#endif
