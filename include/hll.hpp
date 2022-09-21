#ifndef HEXED_HLL_HPP_
#define HEXED_HLL_HPP_

#include <algorithm>
#include <cmath>

namespace hexed
{

/*
 * Computes the numerical face flux using a Harten-Lax-van Leer (HLL) flux.
 * The face state is obtained from `data` (layout: `[i_side][i_var][i_face_qpoint]`).
 * The resulting numerical flux is then written to `data` (with a copy for each side).
 */
template<int n_dim, int n_face_qpoint>
void hll(double* data, int i_dim, double heat_rat)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;

  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double flux[2][n_var];
    double velocity[2];
    double sound_speed[2];

    // compute flux and other variables on each side
    for (int i_side = 0; i_side < 2; ++i_side) {
      #define READ(i) data[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
      #define FLUX(i) flux[i_side][i]
      double veloc = READ(i_dim)/READ(n_var - 2);
      double pres = 0;
      for (int j_dim = 0; j_dim < n_var - 2; ++j_dim) {
        FLUX(j_dim) = READ(j_dim)*veloc;
        pres += READ(j_dim)*READ(j_dim)/READ(n_var - 2);
      }
      pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
      FLUX(i_dim) += pres;
      FLUX(n_var - 2) = READ(i_dim);
      FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
      velocity[i_side] = veloc;
      sound_speed[i_side] = std::sqrt(heat_rat*pres/READ(n_var - 2));
      #undef FLUX
      #undef READ
    }
    double wave_speed [2];
    wave_speed[0] = std::min(velocity[0] - sound_speed[0], velocity[1] - sound_speed[1]);
    wave_speed[1] = std::max(velocity[0] + sound_speed[0], velocity[1] + sound_speed[1]);

    // get numerical flux from the flux/state on both sides
    for (int i_var = 0; i_var < n_var; ++i_var) {
      const int i = i_var*n_face_qpoint + i_qpoint;
      double num_flux;
      if      (wave_speed[0] >= 0) num_flux = flux[0][i_var];
      else if (wave_speed[1] <= 0) num_flux = flux[1][i_var];
      else {
        num_flux = (wave_speed[1]*flux[0][i_var] - wave_speed[0]*flux[1][i_var]
                    + wave_speed[0]*wave_speed[1]*(data[i + face_size] - data[i]))
                   / (wave_speed[1] - wave_speed[0]);
      }
      data[i            ] = num_flux;
      data[i + face_size] = num_flux;
    }
  }
}

}
#endif
