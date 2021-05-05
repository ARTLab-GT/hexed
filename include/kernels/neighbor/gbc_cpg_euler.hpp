#ifndef CARTDG_GBC_CPG_EULER_HPP_
#define CARTDG_GBC_CPG_EULER_HPP_

#include "cpg_euler_hll_deformed.hpp"
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "../Kernel_settings.hpp"
#include "../../Ghost_boundary_condition.hpp"

namespace cartdg
{

// AUTOGENERATE
template<int n_var, int n_qpoint, int row_size>
void gbc_cpg_euler(std::vector<Ghost_boundary_condition*>& ghost_bound_conds,
                   double* read, double* write, double weight,
                   Kernel_settings& settings)
{
  const int n_dof = n_var*n_qpoint;
  const int n_dim = n_var - 2;

  for (Ghost_boundary_condition* gbc : ghost_bound_conds)
  {
    for (unsigned i_bc = 0; i_bc < gbc->elems.size(); ++i_bc)
    {
      int i_elem = gbc->elems[i_bc];
      int stride = n_qpoint;
      for (int i_axis = 0; i_axis < gbc->i_dim + 1; ++i_axis) stride /= row_size;
      read_copy<n_var, n_qpoint, row_size>(read + i_elem*n_dof, gbc->domain_state().data(), stride, gbc->is_positive_face);
      const int jac_size = n_dim*n_dim*n_qpoint/row_size;
      double face_jacobian [2][jac_size];
      read_copy<n_dim*n_dim, n_qpoint, row_size>(gbc->jacobians[i_bc], &face_jacobian[0][0], stride, gbc->is_positive_face);
      for (int i = 0; i < jac_size; ++i)
      {
        face_jacobian[1][i] = face_jacobian[0][i];
      }
      double mult = settings.d_t_by_d_pos/weight;
      gbc->calc_ghost_state();
      Eigen::ArrayXXd d_flux (n_qpoint/row_size, 2*n_var);
      int i_axis_arg [] {gbc->i_dim, gbc->i_dim};
      bool flip [] {false, false};
      cpg_euler_hll_deformed<n_dim, n_qpoint/row_size>(gbc->state.data(), d_flux.data(), &face_jacobian[0][0], mult, i_axis_arg, flip, settings.cpg_heat_rat);
      double* d_flux_r = d_flux.data();
      if (!gbc->is_positive_face) d_flux_r += n_qpoint/row_size*n_var;
      write_copy<n_var, n_qpoint, row_size>(d_flux_r, write + i_elem*n_dof, stride, gbc->is_positive_face);
    }
  }
}

}

#endif
