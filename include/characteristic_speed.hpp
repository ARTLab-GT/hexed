#ifndef CARTDG_CHARACTERISTIC_SPEED_HPP_
#define CARTDG_CHARACTERISTIC_SPEED_HPP_

#include <cmath>

namespace cartdg
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
  return max_speed;
}

}
#endif
