#ifndef CARTDG_GBC_CONVECTIVE_HPP_
#define CARTDG_GBC_CONVECTIVE_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Ghost_boundary_condition.hpp>
#include <Grid.hpp>
#include "hll_deformed_cpg_euler.hpp"
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void gbc_convective(Grid& grid, Basis& basis, Kernel_settings& settings)
{
  double weight = basis.node_weights()[0];
  const int n_dim = n_var - 2;
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;

  for (Ghost_boundary_condition* gbc : grid.ghost_bound_conds)
  {
    //#pragma omp parallel for
    for (unsigned i_bc = 0; i_bc < gbc->elems.size(); ++i_bc)
    {
      int i_elem = gbc->elems[i_bc];
      Element& elem {grid.element(i_elem)};
      int stride = n_qpoint;
      for (int i_dim = 0; i_dim < gbc->i_dim + 1; ++i_dim) stride /= row_size;
      read_copy<n_var, n_qpoint, row_size>(elem.stage(i_read), gbc->domain_state().data(), stride, gbc->is_positive_face);
      const int jac_size = n_dim*n_dim*n_qpoint/row_size;
      double face_jacobian [2][jac_size];
      {
        double jacobian[n_dim][n_dim][n_qpoint];
        for (int i_jac = 0; i_jac < n_dim; ++i_jac) // this is very inefficient. if it becomes a problem it can certainly be optimized
        for (int j_jac = 0; j_jac < n_dim; ++j_jac)
        {
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
          {
            jacobian[i_jac][j_jac][i_qpoint] = elem.jacobian(i_jac, j_jac, i_qpoint);
          }
        }
        read_copy<n_dim*n_dim, n_qpoint, row_size>(&jacobian[0][0][0], &face_jacobian[0][0], stride, gbc->is_positive_face);
      }
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
      write_copy<n_var, n_qpoint, row_size, true>(d_flux_r,elem.stage(i_write), stride, gbc->is_positive_face);
    }
  }
}

}

#endif
