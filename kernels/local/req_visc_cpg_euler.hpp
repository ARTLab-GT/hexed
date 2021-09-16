#ifndef CARTDG_REQ_VISC_CPG_EULER_HPP_
#define CARTDG_REQ_VISC_CPG_EULER_HPP_

#include <math.hpp>
#include <Basis.hpp>
#include <Kernel_settings.hpp>
#include "../observing/indicator.hpp"
#include "../observing/char_speed_convective.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(regular, 1)
template<int n_var, int n_qpoint, int row_size>
void req_visc_cpg_euler(double* read, double* visc, int n_elem, Basis& basis, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  constexpr int n_visc = custom_math::pow(2, n_dim);
  const double heat_rat = settings.cpg_heat_rat;
  const double d_pos = settings.d_pos;
  double weights [row_size];
  double ortho [row_size];
  {
    auto weights_mat = basis.node_weights();
    auto ortho_mat = basis.orthogonal(row_size - 1);
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
    {
      weights[i_qpoint] = weights_mat(i_qpoint);
      ortho[i_qpoint] = ortho_mat(i_qpoint);
    }
  }

  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    double* elem_read = read + i_elem*n_qpoint*n_var;
    double mass_indicator = indicator<n_qpoint, row_size>(elem_read + n_dim*n_qpoint, weights, ortho);
    double char_speed = char_speed_convective<n_var, n_qpoint>(elem_read, heat_rat);
    double req_visc = mass_indicator*char_speed*d_pos/(row_size - 1.);
    for (int i_visc = 0; i_visc < n_visc; ++i_visc)
    {
      visc[i_elem*n_visc + i_visc] = req_visc;
    }
  }
}

}
#endif
