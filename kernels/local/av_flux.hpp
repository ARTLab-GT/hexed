#ifndef CARTDG_AV_FLUX_HPP_
#define CARTDG_AV_FLUX_HPP_

#include <Basis.hpp>
#include <Kernel_settings.hpp>
#include "../static_math.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void av_flux(double* flux, double* visc, int n_elem, Basis& basis, Kernel_settings& settings)
{
  double d_t = settings.d_t_by_d_pos*settings.d_pos;
  const int n_dim = n_var - 2;
  constexpr int n_visc = static_math::pow(2, n_dim);

  double interp [n_visc][n_qpoint];
  for (int i_visc = 0; i_visc < n_visc; ++i_visc)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      interp[i_visc][i_qpoint] = 1.;
    }
  }
  for (int stride_f = n_qpoint/row_size, stride_v = n_visc/2, n_rows = 1, i_axis = 0;
       n_rows < n_qpoint;
       stride_f /= row_size, stride_v /= 2, n_rows *= row_size, ++i_axis)
  {
    for (int i_outer = 0; i_outer < n_rows; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride_f; ++i_inner)
      {
        for (int i_visc = 0; i_visc < n_visc; ++i_visc)
        {
          bool positive = ((i_visc/stride_v)%2) == 1;
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            double node = basis.node(i_qpoint);
            interp[i_visc][i_outer*stride_f*row_size + i_inner + i_qpoint*stride_f]
            *= (positive ? node : 1. - node);
          }
        }
      }
    }
  }

  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      double qpoint_visc = 0.;
      for (int i_visc = 0; i_visc < n_visc; ++i_visc)
      {
        qpoint_visc += visc[i_elem*n_visc + i_visc]*interp[i_visc][i_qpoint];
      }
      flux[i_elem*n_qpoint + i_qpoint] *= qpoint_visc*d_t;
    }
  }
}

}
#endif
