#ifndef CARTDG_MCS_CONVECTIVE_HPP_
#define CARTDG_MCS_CONVECTIVE_HPP_

#include <Vector_view.hpp>
#include <Kernel_settings.hpp>
#include <Element.hpp>
#include "char_speed_convective.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(cartesian, 1)
template<int n_var, int n_qpoint, int row_size>
double mcs_convective(Sequence<Element&>& elements, Kernel_settings& settings)
{
  double heat_rat = settings.cpg_heat_rat;
  double max_speed = 0.;
  #pragma omp parallel for reduction(max:max_speed)
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double max_tss = 0.;
    double* tss = elements[i_elem].time_step_scale();
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      max_tss = std::max(max_tss, tss[i_qpoint]);
    }
    max_speed = std::max(char_speed_convective<n_var, n_qpoint>(elements[i_elem].stage(0), heat_rat)*max_tss, max_speed);
  }
  return max_speed;
}

}
#endif
