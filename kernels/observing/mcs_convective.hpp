#ifndef CARTDG_MCS_CONVECTIVE_HPP_
#define CARTDG_MCS_CONVECTIVE_HPP_

#include <Kernel_settings.hpp>
#include <Element.hpp>
#include "char_speed_convective.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(regular, 1)
template<int n_var, int n_qpoint, int row_size>
double mcs_convective(elem_vec& elements, Kernel_settings& settings)
{
  double heat_rat = settings.cpg_heat_rat;
  const int i_read = settings.i_read;
  double max_speed = 0.;
  #pragma omp parallel for reduction(max:max_speed)
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    max_speed = std::max<double>(char_speed_convective<n_var, n_qpoint>(elements[i_elem]->stage(i_read), heat_rat), max_speed);
  }
  return max_speed;
}

}
#endif
