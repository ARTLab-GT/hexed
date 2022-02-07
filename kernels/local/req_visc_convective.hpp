#ifndef CARTDG_REQ_VISC_CONVECTIVE_HPP_
#define CARTDG_REQ_VISC_CONVECTIVE_HPP_

#include <memory>

#include <math.hpp>
#include <Basis.hpp>
#include <Kernel_settings.hpp>
#include <Deformed_element.hpp>
#include "../observing/indicator.hpp"
#include "../observing/char_speed_convective.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size, typename E>
void req_visc_convective(std::vector<std::unique_ptr<E>>& elements, Basis& basis, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  constexpr int n_visc = custom_math::pow(2, n_dim);
  const double heat_rat = settings.cpg_heat_rat;
  const double d_pos = settings.d_pos;
  const int i_read = settings.i_read;
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

  for (std::unique_ptr<E>& elem : elements)
  {
    double* stage = elem->stage(i_read);
    double mass_indicator = indicator<n_qpoint, row_size>(stage + n_dim*n_qpoint, weights, ortho);
    double char_speed = char_speed_convective<n_var, n_qpoint>(stage, elem->time_step_scale(), heat_rat);
    double req_visc = mass_indicator*char_speed*d_pos/(row_size - 1.);
    double* visc = elem->viscosity();
    for (int i_visc = 0; i_visc < n_visc; ++i_visc)
    {
      visc[i_visc] = req_visc;
    }
  }
}

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void req_visc_regular_convective(elem_vec& elements, Basis& basis, Kernel_settings& settings)
{
  req_visc_convective<n_var, n_qpoint, row_size, Element>(elements, basis, settings);
}

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void req_visc_deformed_convective(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  req_visc_convective<n_var, n_qpoint, row_size, Deformed_element>(def_elements, basis, settings);
}

}
#endif
