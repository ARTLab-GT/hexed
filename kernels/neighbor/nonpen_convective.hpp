#ifndef CARTDG_NONPEN_CONVECTIVE_HPP_
#define CARTDG_NONPEN_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template <int n_var, int n_qpoint, int row_size>
void nonpen_convective(def_elem_wall_vec& walls, Basis& basis, Kernel_settings& settings)
{
  const double heat_rat = settings.cpg_heat_rat;
  const int n_dim = n_var - 2;
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_var*n_face_qpoint;

  #pragma omp parallel for
  for (unsigned i_bc = 0; i_bc < walls.size(); ++i_bc)
  {
    auto face_index {walls[i_bc].face_index()};
    const int i_dim = face_index.i_dim;
    int is_p = face_index.is_positive;
    Element& elem {*face_index.element};

    double both_face [2][n_var][n_face_qpoint];
    double* dom_face = elem.face() + (i_dim*2 + is_p)*face_size;
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
      {
        both_face[0][i_var][i_qpoint] = both_face[1][i_var][i_qpoint] = dom_face[i_var*n_face_qpoint + i_qpoint];
      }
    }
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
    {
      both_face[is_p][i_dim][i_qpoint] *= -1;
    }
    hll_cpg_euler<n_dim, n_face_qpoint>(both_face[0][0], both_face[0][0], 1., i_dim, heat_rat);
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
      {
        dom_face[i_var*n_face_qpoint + i_qpoint] = both_face[1 - is_p][i_var][i_qpoint];
      }
    }
  }
}

}

#endif
