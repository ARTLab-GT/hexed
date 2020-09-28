#ifndef CPG_EULER_HLL_HPP_
#define CPG_EULER_HLL_HPP_

template<int n_dim, int n_face_qpoint>
void cpg_euler_hll(double* state_r, double* d_flux_w, double mult,
                   int i_axis, double sp_heat_rat=1.4)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;
  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double flux[n_var*2];
    double wave_speed[2];
    for (int i_side = 0; i_side < 2; ++i_side)
    {
      #define READ(i) state_r[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
      #define FLUX(i) flux[i + n_var*i_side]
      double veloc = READ(i_axis)/READ(n_var - 2);
      double pres = 0;
      for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
      {
        FLUX(j_axis) = READ(j_axis)*veloc;
        pres += READ(j_axis)*READ(j_axis)/READ(n_var - 2);
      }
      pres = (sp_heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
      FLUX(i_axis) += pres;
      FLUX(n_var - 2) = READ(i_axis);
      FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
      wave_speed[i_side] = veloc + (2*i_side - 1)*sqrt(sp_heat_rat*pres/READ(n_var - 2));
      #undef FLUX
      #undef READ
    }
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      double num_flux;
      if      (wave_speed[0] >= 0) num_flux = flux[i_var];
      else if (wave_speed[1] <= 0) num_flux = flux[i_var + n_var];
      else
      {
        num_flux = (  wave_speed[1]*flux[i_var        ]
                    - wave_speed[0]*flux[i_var + n_var]
                    + wave_speed[0]*wave_speed[1]*(
                      state_r[i_var*n_face_qpoint + i_qpoint + face_size]
                    - state_r[i_var*n_face_qpoint + i_qpoint]))
                   / (wave_speed[1] - wave_speed[0]);
      }
      d_flux_w[i_qpoint + i_var*n_face_qpoint            ]
      = (flux[i_var] - num_flux)*mult;
      d_flux_w[i_qpoint + i_var*n_face_qpoint + face_size]
      = (num_flux - flux[i_var + n_var])*mult;
    }
  }
}

#endif
