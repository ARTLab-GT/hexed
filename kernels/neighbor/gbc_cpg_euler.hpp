#ifndef CARTDG_GBC_CPG_EULER_HPP_
#define CARTDG_GBC_CPG_EULER_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Ghost_boundary_condition.hpp>
#include "hll_deformed_cpg_euler.hpp"
#include "read_copy.hpp"
#define CARTDG_ATOMIC
#include "write_copy.hpp"
#undef CARTDG_ATOMIC

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void gbc_cpg_euler(std::vector<Ghost_boundary_condition*>& ghost_bound_conds,
                   double* read, double* write, Basis& basis,
                   Kernel_settings& settings)
{
  double weight = basis.node_weights()[0];
  const int n_dof = n_var*n_qpoint;
  const int n_dim = n_var - 2;

  #pragma omp parallel for
  for (unsigned i_gbc = 0; i_gbc < ghost_bound_conds.size(); ++i_gbc) // range-based for not supported by omp <=4.0
  {
    Ghost_boundary_condition* gbc = ghost_bound_conds[i_gbc];
    for (unsigned i_bc = 0; i_bc < gbc->elems.size(); ++i_bc)
    {
      int i_elem = gbc->elems[i_bc];
      int stride = n_qpoint;
      for (int i_dim = 0; i_dim < gbc->i_dim + 1; ++i_dim) stride /= row_size;
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
      int i_dim_arg [] {gbc->i_dim, gbc->i_dim};
      bool flip [] {false, false};
      hll_deformed_cpg_euler<n_dim, n_qpoint/row_size>(gbc->state.data(), d_flux.data(), &face_jacobian[0][0], mult, i_dim_arg, flip, settings.cpg_heat_rat);
      double* d_flux_r = d_flux.data();
      if (!gbc->is_positive_face) d_flux_r += n_qpoint/row_size*n_var;
      write_copy<n_var, n_qpoint, row_size>(d_flux_r, write + i_elem*n_dof, stride, gbc->is_positive_face);
    }
  }
}

}

#endif
