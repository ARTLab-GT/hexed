#ifndef CPG_EULER_FBC_HPP_
#define CPG_EULER_FBC_HPP_

#include "cpg_euler_hll.hpp"
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "../Kernel_settings.hpp"
#include "../../Fitted_boundary_condition.hpp"

namespace cartdg
{

template <int n_var, int n_qpoint, int row_size>
void cpg_euler_fbc(std::vector<Fitted_boundary_condition*>& fit_bound_conds,
                   double* read, double* write,
                   double weight,
                   Kernel_settings& settings)
{
  const int n_dof = n_var*n_qpoint;

  for (Fitted_boundary_condition* fbc : fit_bound_conds)
  {
    for (int i_elem : fbc->elems)
    {
      int stride = n_qpoint;
      for (int i_axis = 0; i_axis < fbc->i_dim + 1; ++i_axis) stride /= row_size;
      read_copy<n_var, n_qpoint, row_size>(read + i_elem*n_dof, fbc->domain_state().data(),
                                           stride, fbc->is_positive_face);
      double mult = settings.d_t_by_d_pos/weight;
      fbc->calc_ghost_state();
      // FIXME: get correct specific heat ratio
      Eigen::ArrayXXd d_flux (n_qpoint/row_size, 2*n_var);
      cpg_euler_hll<n_var - 2, n_qpoint/row_size>(fbc->state.data(), d_flux.data(), mult,
                                                  fbc->i_dim, settings.cpg_heat_rat);
      double* d_flux_r = d_flux.data();
      if (!fbc->is_positive_face) d_flux_r += n_qpoint/row_size*n_var;
      write_copy<n_var, n_qpoint, row_size>(d_flux_r, write + i_elem*n_dof,
                                            stride, fbc->is_positive_face);
    }
  }
}

}

#endif