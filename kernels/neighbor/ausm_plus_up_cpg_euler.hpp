#ifndef CARTDG_AUSM_PLUS_UP_CPG_EULER
#define CARTDG_AUSM_PLUS_UP_CPG_EULER

#include <limits>

namespace cartdg
{

const double ref_mach = 0.3;

template<int n_dim, int n_face_qpoint>
void ausm_plus_up_cpg_euler(double* state_r, double* d_flux_w, double mult,
                            int i_dim, double heat_rat)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;
  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double flux [2][n_var];
    double mid_sound_speed = std::numeric_limits<double>::max();
    double veloc [2];
    double pres [2];
    double mass [2];
    for (int i_side : {0, 1})
    {
      #define READ(i) state_r[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
      #define FLUX(i) flux[i_side][i]
      mass[i_side] = READ(n_var - 2);
      veloc[i_side] = READ(i_dim)/READ(n_var - 2);
      pres[i_side] = 0;
      for (int j_dim = 0; j_dim < n_var - 2; ++j_dim)
      {
        FLUX(j_dim) = READ(j_dim)*veloc[i_side];
        pres[i_side] += READ(j_dim)*READ(j_dim)/READ(n_var - 2);
      }
      pres[i_side] = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres[i_side]);
      FLUX(i_dim) += pres[i_side];
      FLUX(n_var - 2) = READ(i_dim);
      FLUX(n_var - 1) = (READ(n_var - 1) + pres[i_side])*veloc[i_side];
      #undef FLUX

      double sonic_sound_speed = std::sqrt(2*(heat_rat - 1)/(heat_rat + 1)*(READ(n_var - 1) + pres[i_side])/READ(n_var - 2));
      mid_sound_speed = std::min(sonic_sound_speed*sonic_sound_speed/std::max(sonic_sound_speed, veloc[i_side]*(1 - 2*i_side)), mid_sound_speed);
      #undef READ
    }
    double mid_mass = 0.5*(mass[0] + mass[1]);
    double mean_sq_mach = (veloc[0]*veloc[0] + veloc[1]*veloc[1])/(2*mid_sound_speed*mid_sound_speed);
    double bounded_mach = std::sqrt(std::min(1., std::max(mean_sq_mach, ref_mach*ref_mach)));
    double sound_mult = bounded_mach*(2. - bounded_mach);

    double mach_poly1 [2][2];
    double mach_poly2 [2][2];
    double mach_poly4 [2][2];
    double pres_poly5 [2][2];
    for (int i_side : {0, 1})
    {
      #define READ(i) state_r[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
      double mach = veloc[i_side]/mid_sound_speed;
      for (int i_sign : {0, 1})
      {
        int sign = 2*i_sign - 1;
        mach_poly1[i_side][i_sign] = 0.5*(mach + sign*std::abs(mach));
        mach_poly2[i_side][i_sign] = sign*0.25*(mach + sign)*(mach + sign);
      }
      for (int i_sign : {0, 1})
      {
        int sign = 2*i_sign - 1;
        mach_poly4[i_side][i_sign] = (std::abs(mach) >= 1) ? mach_poly1[i_side][i_sign] : mach_poly2[i_side][i_sign]*(1. - sign*2*mach_poly2[i_side][1 - i_sign]);
        pres_poly5[i_side][i_sign] = (std::abs(mach) >= 1) ? mach_poly1[i_side][i_sign]/mach : mach_poly2[i_side][i_sign]*(sign*2 - mach - sign*3*(-4. + 5*sound_mult*sound_mult)*mach*mach_poly2[i_side][1 - i_sign]);
      }
      #undef READ
    }
    double mid_mach = mach_poly4[0][1] + mach_poly4[1][0] - 0.25/sound_mult*std::max(1. - mean_sq_mach, 0.)*(pres[1] - pres[0])/(mid_mass*mid_sound_speed*mid_sound_speed);
    double mass_flux = mid_sound_speed*mid_mach*((mid_mach > 0) ? mass[0] : mass[1]);
    double mid_pres = pres_poly5[0][1]*pres[0] + pres_poly5[1][0]*pres[1] - 0.75*pres_poly5[0][1]*pres_poly5[1][0]*(mass[0] + mass[1])*sound_mult*mid_sound_speed*(veloc[1] - veloc[0]);

    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      const int i = i_var*n_face_qpoint + i_qpoint;
      int upwind_ind = mass_flux > 0 ? 0 : 1;
      int upwind_offset = upwind_ind*face_size;
      double num_flux = state_r[i + upwind_offset];
      if (i_var == n_var - 1) num_flux += pres[upwind_ind];
      num_flux *= mass_flux/mass[upwind_ind];
      if (i_var == i_dim) num_flux += mid_pres;
      d_flux_w[i            ] = num_flux;
      d_flux_w[i + face_size] = num_flux;
    }
  }
}

}

#endif
