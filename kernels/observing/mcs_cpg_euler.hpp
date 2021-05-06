#ifndef CARTDG_MCS_CPG_EULER_HPP_
#define CARTDG_MCS_CPG_EULER_HPP_

#include <cmath>
#include <stdexcept>

#include <Kernel_settings.hpp>

namespace cartdg
{

// AUTOGENERATE
template<int n_var, int n_qpoint, int row_size>
double mcs_cpg_euler(double* read, int n_elem, Kernel_settings& settings)
{
  const int n_dof = n_var*n_qpoint;
  const int n_dim = n_var - 2;
  double heat_rat = settings.cpg_heat_rat;

  double max_speed = 0.;
  #pragma omp parallel for reduction(max:max_speed)
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      #define READ(i) read[n_dof*i_elem + (i)*n_qpoint + i_qpoint]
      double max_veloc = 0;
      double int_ener = 0;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        const double veloc = READ(i_dim)/READ(n_var - 2);
        max_veloc = std::max<double>(max_veloc, abs(veloc));
        int_ener += veloc*READ(i_dim);
      }
      int_ener = READ(n_var - 1) - 0.5*int_ener;
      const double sound_speed = std::sqrt(heat_rat*(heat_rat - 1)
                                           *int_ener/READ(n_var - 2));
      max_speed = std::max(max_speed, max_veloc + sound_speed);
      #undef READ
    }
  }
  return max_speed;
}

}
#endif
