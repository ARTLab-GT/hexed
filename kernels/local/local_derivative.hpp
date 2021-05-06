#ifndef CARTDG_LOCAL_DERIVATIVE_HPP_
#define CARTDG_LOCAL_DERIVATIVE_HPP_

#include "derivative.hpp"

namespace cartdg
{

// AUTOGENERATE
template<int n_var, int n_qpoint, int row_size>
void local_derivative(double* read, double* write, int n_elem, int i_var, int i_axis,
                      Basis& basis, Kernel_settings& settings)
{
  derivative<n_var, 1, n_qpoint, row_size, false>(read, write, n_elem, i_var, 0, i_axis, basis, settings);
}

}
#endif
