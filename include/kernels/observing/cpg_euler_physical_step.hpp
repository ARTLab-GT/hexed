#ifndef CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_
#define CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_

#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var>
double cpg_euler_physical_step(double* read, double* write, int n_points,
                               Kernel_settings& settings)
{
  return 1.;
}

}

#endif