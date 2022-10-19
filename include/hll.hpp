#ifndef HEXED_HLL_HPP_
#define HEXED_HLL_HPP_

#include "thermo.hpp"

namespace hexed::hll
{

typedef std::array<double, 2> arr2;

// core formula of Harten-Lax-van Leer (HLL) flux.
inline double hll(arr2 speed, arr2 flux, arr2 state)
{
  if (speed[0] >= 0) return flux[0];
  if (speed[1] <= 0) return flux[1];
  return (speed[1]*flux[0] - speed[0]*flux[1] + speed[0]*speed[1]*(state[1] - state[0]))
         /(speed[1] - speed[0]);
}

/*
 * Computes the numerical face flux for inviscid flow.
 * The face state is obtained from `data` (layout: `[i_side][i_var][i_face_qpoint]`).
 * The resulting numerical flux is then written to `data` (with a copy for each side).
 */
template<int n_dim, int n_face_qpoint>
void inviscid(double* data, const double* normal, double heat_rat)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;

  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double flux[2][n_var];
    arr2 velocity;
    arr2 sound_speed;
    double nrml [3];
    double nrml_mag = 0.;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      nrml[i_dim] = normal[i_dim*n_face_qpoint + i_qpoint];
      nrml_mag += nrml[i_dim]*nrml[i_dim];
    }
    nrml_mag = std::sqrt(nrml_mag);

    // compute flux and other variables on each side
    // note: flux, velocity, sound speed, and wave speed are in **reference space**
    for (int i_side = 0; i_side < 2; ++i_side) {
      #define READ(i) data[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
      #define FLUX(i) flux[i_side][i]
      HEXED_COMPUTE_SCALARS
      HEXED_ASSERT_ADMISSIBLE
      double veloc = 0;
      for (int j_dim = 0; j_dim < n_var - 2; ++j_dim) {
        veloc += READ(j_dim)/mass*nrml[j_dim];
      }
      for (int j_dim = 0; j_dim < n_var - 2; ++j_dim) {
        FLUX(j_dim) = READ(j_dim)*veloc + pres*nrml[j_dim];
      }
      FLUX(n_var - 2) = mass*veloc;
      FLUX(n_var - 1) = (ener + pres)*veloc;
      velocity[i_side] = veloc;
      sound_speed[i_side] = std::sqrt(heat_rat*pres/mass)*nrml_mag;
      HEXED_ASSERT(!std::isnan(sound_speed[i_side]), "speed of sound is NaN")
      #undef FLUX
      #undef READ
    }
    arr2 wave_speed;
    wave_speed[0] = std::min(velocity[0] - sound_speed[0], velocity[1] - sound_speed[1]);
    wave_speed[1] = std::max(velocity[0] + sound_speed[0], velocity[1] + sound_speed[1]);

    // get numerical flux from the flux/state on both sides
    for (int i_var = 0; i_var < n_var; ++i_var) {
      const int i = i_var*n_face_qpoint + i_qpoint;
      data[i] = data[i + face_size] = hll(wave_speed, {flux[0][i_var], flux[1][i_var]}, {data[i], data[i + face_size]});
    }
  }
}

/*
 * Computes the numerical face flux for linear advection.
 * The first `n_dim` variables of `data` contain the advection velocity.
 * Variable `n_dim` contains the scalar state.
 * The scalar flux is written to variable `n_dim` on both sides.
 * The velocity is left unchanged.
 */
template<int n_dim, int n_face_qpoint>
void advection(double* data, const double* normal)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;

  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double nrml [3];
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      nrml[i_dim] = normal[i_dim*n_face_qpoint + i_qpoint];
    }
    // compute flux on each side
    // note: flux and velocity are in **reference space**
    arr2 veloc {0., 0.};
    arr2 flux;
    arr2 state;
    for (int i_side = 0; i_side < 2; ++i_side) {
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        veloc[i_side] += data[(i_dim + i_side*n_var)*n_face_qpoint + i_qpoint]*nrml[i_dim];
      }
      state[i_side] = data[(n_dim + i_side*n_var)*n_face_qpoint + i_qpoint];
      flux[i_side] = state[i_side]*veloc[i_side];
    }
    arr2 wave_speed;
    for (int i_side = 0; i_side < 2; ++i_side) {
      int sign = 2*i_side - 1;
      wave_speed[i_side] = std::min(veloc[0] + sign*.01, veloc[1] + sign*.01);
    }
    // get numerical flux from the flux/state on both sides
    double* scalar_flux = data + 2*n_face_qpoint + i_qpoint;
    scalar_flux[0] = scalar_flux[face_size] = hll(wave_speed, flux, state);
  }
}

}
#endif
