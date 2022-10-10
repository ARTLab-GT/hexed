#ifndef HEXED_CHARACTERISTIC_SPEED_HPP_
#define HEXED_CHARACTERISTIC_SPEED_HPP_

#include <cmath>
#include "thermo.hpp"

namespace hexed
{

/*
 * Computes the maximum speed of characteristics in physical space,
 * not accounting for local time stepping or any other numerical artifiacts.
 */
template<int n_dim, int n_qpoint>
double characteristic_speed(double* read, double heat_rat)
{
  const int n_var = n_dim + 2;
  double max_speed = 0.;
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
  {
    #define READ(i) read[(i)*n_qpoint + i_qpoint]
    HEXED_COMPUTE_SCALARS
    HEXED_ASSERT_ADMISSIBLE
    const double sound_speed = std::sqrt(heat_rat*pres/mass);
    const double veloc = std::sqrt(mmtm_sq/mass/mass);
    max_speed = std::max(max_speed, sound_speed + veloc);
    #undef READ
  }
  return max_speed;
}

}
#endif
