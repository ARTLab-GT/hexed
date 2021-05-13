#ifndef CARTDG_MCS_CPG_EULER_HPP_
#define CARTDG_MCS_CPG_EULER_HPP_

#include <Kernel_settings.hpp>
#include "char_speed_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK
template<int n_var, int n_qpoint, int row_size>
double mcs_cpg_euler(double* read, int n_elem, Kernel_settings& settings)
{
  const int n_dof = n_var*n_qpoint;
  double heat_rat = settings.cpg_heat_rat;
  double max_speed = 0.;
  #pragma omp parallel for reduction(max:max_speed)
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    max_speed = std::max<double>(char_speed_cpg_euler<n_var, n_qpoint>(read + n_dof*i_elem, heat_rat), max_speed);
  }
  return max_speed;
}

}
#endif
