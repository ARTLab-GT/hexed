#ifndef CARTDG_CPG_EULER_RESTRICT_STEP_HPP_
#define CARTDG_CPG_EULER_RESTRICT_STEP_HPP_

#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_restrict_step(double* read, double* write, int n_elem,
                             double step, Kernel_settings& settings)
{
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i = 0; i < n_var*n_qpoint; ++i)
    {
      double& w = write[n_var*n_qpoint*i_elem + i];
      double& r = read[n_var*n_qpoint*i_elem + i];
      w = r + step*(w - r);
    }
  }
}

}

#endif