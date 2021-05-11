#ifndef CARTDG_GBC_AV_HPP_
#define CARTDG_GBC_AV_HPP_

#include <Kernel_settings.hpp>
#include <Ghost_boundary_condition.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template <int n_var, int n_qpoint, int row_size>
void gbc_av(std::vector<Ghost_boundary_condition*>& ghost_bound_conds,
            double* read, double* write,
            int i_var, int i_axis, double weight,
            Kernel_settings& settings)
{
  double mult = 1./(weight*settings.d_pos);
  const int n_dof = n_var*n_qpoint;
  for (Ghost_boundary_condition* gbc : ghost_bound_conds)
  {
    if (gbc->i_dim == i_axis)
    {
      int stride = n_qpoint;
      for (int i_axis = 0; i_axis < gbc->i_dim + 1; ++i_axis) stride /= row_size;
      for (unsigned i_bc = 0; i_bc < gbc->elems.size(); ++i_bc)
      {
        int i_elem = gbc->elems[i_bc];
        double flux [n_qpoint] {};
        read_copy<1, n_qpoint, row_size>(read + i_elem*n_qpoint, flux, stride, gbc->is_positive_face);
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
        {
          flux[i_qpoint] *= mult;
          if (gbc->is_positive_face) flux[i_qpoint] *= -1;
        }
        write_copy<1, n_qpoint, row_size>(flux, write + i_elem*n_dof + i_var*n_qpoint, stride, gbc->is_positive_face);
      }
    }
  }
}

}

#endif
