#ifndef CARTDG_GBC_AV_HPP_
#define CARTDG_GBC_AV_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Ghost_boundary_condition.hpp>
#include <Grid.hpp>
#include <math.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template <int n_var, int n_qpoint, int row_size>
void gbc_av(Grid& grid, int i_var, int i_dim, Basis& basis, Kernel_settings& settings)
{
  #if 0
  const int n_dim = n_var - 2;
  const int i_write = settings.i_write;
  const int stride = custom_math::pow(row_size, n_dim - i_dim - 1);
  double mult = 1./(basis.node_weights()[0]*settings.d_pos);

  for (Ghost_boundary_condition* gbc : grid.ghost_bound_conds)
  {
    if (gbc->i_dim == i_dim)
    {
      //#pragma omp parallel for
      for (unsigned i_bc = 0; i_bc < gbc->elems.size(); ++i_bc)
      {
        int i_elem = gbc->elems[i_bc];
        Element& elem {grid.element(i_elem)};
        double flux [n_qpoint/row_size] {};
        read_copy<1, n_qpoint, row_size>(elem.derivative(), flux, stride, gbc->is_positive_face);
        for (int i_qpoint = 0; i_qpoint < n_qpoint/row_size; ++i_qpoint)
        {
          flux[i_qpoint] *= mult;
          if (gbc->is_positive_face) flux[i_qpoint] *= -1;
        }
        write_copy<1, n_qpoint, row_size>(flux, elem.stage(i_write) + i_var*n_qpoint, stride, gbc->is_positive_face);
      }
    }
  }
  #endif
}

}

#endif
