#ifndef CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_
#define CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_

#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint>
double cpg_euler_physical_step(double* read, double* write, int n_elem,
                               Kernel_settings& settings)
{
  double max_diff = settings.max_difference;
  double step = 1.;
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      #define READ(i_var) read[i_elem*n_var*n_qpoint + (i_var)*n_qpoint + i_qpoint]
      #define WRITE(i_var) write[i_elem*n_var*n_qpoint + (i_var)*n_qpoint + i_qpoint]
      double mass_diff = 1. - WRITE(n_var - 2)/READ(n_var - 2);
      if (mass_diff > 0) step = std::min(step, max_diff/mass_diff);
      #undef READ
      #undef WRITE
    }
  }
  return step;
}

}

#endif