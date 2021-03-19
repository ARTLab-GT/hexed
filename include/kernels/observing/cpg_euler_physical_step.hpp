#ifndef CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_
#define CARTDG_CPG_EULER_PHYSICAL_STEP_HPP_

#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
double cpg_euler_physical_step(double* read, double* write, int n_elem,
                               Kernel_settings& settings)
{
  double max_diff = settings.max_difference;
  double step = 1.;
  #pragma omp parallel for reduction(min:step)
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      #define READ(i_var) read[i_elem*n_var*n_qpoint + (i_var)*n_qpoint + i_qpoint]
      #define WRITE(i_var) write[i_elem*n_var*n_qpoint + (i_var)*n_qpoint + i_qpoint]
      double mass_rel_diff = 1. - WRITE(n_var - 2)/READ(n_var - 2);
      if (mass_rel_diff > 0) step = std::min(step, max_diff/mass_rel_diff);

      double ener_diff = WRITE(n_var - 1) - READ(n_var - 1);
      double kin_ener_diff = 0.;
      double kin_ener = 0.;
      for (int i_axis = 0; i_axis < n_var - 2; ++i_axis)
      {
        double diff = 2*(WRITE(i_axis) - READ(i_axis));
        kin_ener_diff += std::max(diff*WRITE(i_axis), diff*READ(i_axis));
        kin_ener += READ(i_axis)*READ(i_axis);
      }
      kin_ener *= 0.5/READ(n_var - 2);
      kin_ener_diff *= 0.5;
      kin_ener_diff = std::max(kin_ener_diff/READ(n_var - 2), kin_ener_diff/WRITE(n_var - 2));
      double pressure_diff = 0.4*(ener_diff - kin_ener_diff);
      double pressure = 0.4*(READ(n_var - 1) - kin_ener);
      double pressure_rel_diff = -pressure_diff/pressure;
      if (pressure_rel_diff > 0) step = std::min(step, max_diff/pressure_rel_diff);
      #undef READ
      #undef WRITE
    }
  }
  return step;
}

}

#endif