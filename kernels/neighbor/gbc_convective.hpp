#ifndef CARTDG_GBC_CONVECTIVE_HPP_
#define CARTDG_GBC_CONVECTIVE_HPP_

#include <iostream>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Ghost_boundary_condition.hpp>
#include <Grid.hpp>
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void gbc_convective(Grid& grid, Basis& basis, Kernel_settings& settings)
{
  const double heat_rat = settings.cpg_heat_rat;
  const int n_dim = n_var - 2;
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_var*n_face_qpoint;

  #pragma omp parallel for // not optimal, but better than nothing
  for (unsigned i_gbc = 0; i_gbc < grid.ghost_bound_conds.size(); ++i_gbc)
  {
    Ghost_boundary_condition* gbc = grid.ghost_bound_conds[i_gbc];
    for (unsigned i_bc = 0; i_bc < gbc->elems.size(); ++i_bc)
    {
      int i_elem = gbc->elems[i_bc];
      Element& elem {grid.element(i_elem)};
      double* face = elem.face() + (gbc->i_dim*2 + gbc->is_positive_face)*face_size;
      double* dom_face = gbc->domain_state().data();
      for (int i_data = 0; i_data < face_size; ++i_data) dom_face[i_data] = face[i_data];
      gbc->calc_ghost_state();
      hll_cpg_euler<n_dim, n_face_qpoint>(gbc->state.data(), gbc->state.data(), 1., gbc->i_dim, heat_rat);
      double* gst_face = gbc->ghost_state().data();
      for (int i_data = 0; i_data < face_size; ++i_data) face[i_data] = gst_face[i_data];
    }
  }
}

}

#endif
