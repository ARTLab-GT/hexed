#ifndef CARTDG_LOCAL_AV_HPP_
#define CARTDG_LOCAL_AV_HPP_

#include "derivative.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_av(double* read, double* write, int n_elem, int i_var, int i_axis,
              Basis& basis, Kernel_settings& settings)
{
  derivative<1, n_var, n_qpoint, row_size, true>(read, write, n_elem, 0, i_var, i_axis, basis, settings);
}

}
#endif
