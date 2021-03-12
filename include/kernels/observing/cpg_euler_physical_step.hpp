#ifndef CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_
#define CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_

#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint>
double cpg_euler_physical_step(double* read, double* write, int n_elem,
                               Kernel_settings& settings)
{
  double step = 1.;
  return step;
}

}

#endif