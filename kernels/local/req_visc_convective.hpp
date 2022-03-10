#ifndef CARTDG_REQ_VISC_CONVECTIVE_HPP_
#define CARTDG_REQ_VISC_CONVECTIVE_HPP_

#include <memory>

#include <math.hpp>
#include <Basis.hpp>
#include <Kernel_settings.hpp>
#include <Deformed_element.hpp>
#include "../observing/indicator.hpp"

namespace cartdg
{

/*
 * Find the artificial viscosity coefficient required in each element and write it to
 * the viscosity at each element vertex. Return the maximum value of the logarithmic
 * nonsmoothness metric to inform smoother convergence criteria.
 */
template<int n_var, int n_qpoint, int row_size, typename E>
double req_visc_convective(std::vector<std::unique_ptr<E>>& elements, Basis& basis, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  constexpr int n_visc = custom_math::pow(2, n_dim);
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

  double max_nonsmooth = -std::numeric_limits<double>::max();
  for (std::unique_ptr<E>& elem : elements)
  {
    double* stage = elem->stage(0);
    auto params = indicator<n_qpoint, row_size>(stage + n_dim*n_qpoint, weights, ortho);
    double mass_indicator = params[0];
    max_nonsmooth = std::max(max_nonsmooth, params[1]);
    double* visc = elem->viscosity();
    for (int i_visc = 0; i_visc < n_visc; ++i_visc)
    {
      visc[i_visc] = mass_indicator;
    }
  }
  return max_nonsmooth;
}

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
double req_visc_regular_convective(elem_vec& elements, Basis& basis, Kernel_settings& settings)
{
  return req_visc_convective<n_var, n_qpoint, row_size, Element>(elements, basis, settings);
}

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
double  req_visc_deformed_convective(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  return req_visc_convective<n_var, n_qpoint, row_size, Deformed_element>(def_elements, basis, settings);
}

}
#endif
