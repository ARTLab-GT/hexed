#ifndef CARTDG_HLL_CPG_EULER_HPP_
#define CARTDG_HLL_CPG_EULER_HPP_

#include <algorithm>
#include <cmath>

namespace cartdg
{

template<int n_dim, int n_face_qpoint>
void hll_cpg_euler(double* state_r, double* d_flux_w, double mult,
                   int i_dim, double sp_heat_rat)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;
  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double flux[2][n_var];
    double velocity[2];
    double sound_speed[2];
    for (int i_side = 0; i_side < 2; ++i_side)
    {
      #define READ(i) state_r[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
      #define FLUX(i) flux[i_side][i]
      double veloc = READ(i_dim)/READ(n_var - 2);
      double pres = 0;
      for (int j_dim = 0; j_dim < n_var - 2; ++j_dim)
      {
        FLUX(j_dim) = READ(j_dim)*veloc;
        pres += READ(j_dim)*READ(j_dim)/READ(n_var - 2);
      }
      pres = (sp_heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
      FLUX(i_dim) += pres;
      FLUX(n_var - 2) = READ(i_dim);
      FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
      velocity[i_side] = veloc;
      sound_speed[i_side] = std::sqrt(sp_heat_rat*pres/READ(n_var - 2));
      #undef FLUX
      #undef READ
    }

    double wave_speed [2];
    wave_speed[0] = std::min(velocity[0] - sound_speed[0], velocity[1] - sound_speed[1]);
    wave_speed[1] = std::max(velocity[0] + sound_speed[0], velocity[1] + sound_speed[1]);

    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      const int i = i_var*n_face_qpoint + i_qpoint;
      double num_flux;
      if      (wave_speed[0] >= 0) num_flux = flux[0][i_var];
      else if (wave_speed[1] <= 0) num_flux = flux[1][i_var];
      else
      {
        num_flux = (wave_speed[1]*flux[0][i_var] - wave_speed[0]*flux[1][i_var]
                    + wave_speed[0]*wave_speed[1]*(state_r[i + face_size] - state_r[i]))
                   / (wave_speed[1] - wave_speed[0]);
      }

      d_flux_w[i            ] = (flux[0][i_var] - num_flux)*mult;
      d_flux_w[i + face_size] = (num_flux - flux[1][i_var])*mult;
    }
  }
}

}
#endif
