#ifndef HEXED_PDE_HPP_
#define HEXED_PDE_HPP_

#include "math.hpp"

/*
 * This namespace contains classes representing the different PDEs Hexed can solve.
 * They are all possible arguments to the `Spatial` class.
 * Each PDE class must have:
 * - `static constexpr int n_var`: the number of state variables required to compute the flux
 *   (some of which may be conserved variables and some of which can be pasive data, e.g. advection velocity)
 * - `static constexpr int curr_start`: index of the variable representing the current state
 * - `static constexpr int ref_start`: index of the variable representing the reference state (for multistage time integration)
 * - `static constexpr int n_update`: the number of conserved variables to be updated
 * - a `flux` function which computes the flux of conserved
 *   variables given a state vector and a (non-unit) reference level normal vector
 * - a `flux_num` function which computes the shared numerical flux given the
 *   state for two neighboring elements and the shared face normal
 */
namespace hexed::pde
{

// computes HLL (Harten-Lax-Van Leer) numerical flux based on wave speed estimate
template <int n_var>
Mat<n_var> hll(Mat<2> speed, Mat<n_var, 2> flux, Mat<n_var, 2> state)
{
  Mat<n_var> f;
  if (speed(0) >= 0) return flux(Eigen::all, 0);
  if (speed(1) <= 0) return flux(Eigen::all, 1);
  return (speed(1)*flux(Eigen::all, 0) - speed(0)*flux(Eigen::all, 1)
          + speed(0)*speed(1)*(state(Eigen::all, 1) - state(Eigen::all, 0)))
         /(speed(1) - speed(0));
}

// contains a PDE class representing the Naver-Stokes equations
// with template options to specify the details of the equation set
template <bool artificial_visc = false>
class Navier_stokes
{
  #define HEXED_ASSERT_THERM_ADMIS \
    HEXED_ASSERT(state(n_dim) > 0, "nonpositive density"); \
    HEXED_ASSERT(state(n_dim + 1) >= 0, "negative energy"); \
    HEXED_ASSERT(pres >= 0, "negative pressure"); \

  public:
  Navier_stokes() = delete;

  template <int n_dim>
  class Pde
  {
    public:
    Pde() = delete;
    static constexpr bool is_viscous = artificial_visc;
    static constexpr bool has_convection = true;
    static constexpr int n_var = n_dim + 2;
    static constexpr int curr_start = 0;
    static constexpr int ref_start = n_var;
    static constexpr int visc_start = 2*n_var;
    static constexpr int n_update = n_var;
    static constexpr double heat_rat = 1.4;
    static constexpr int tss_pow = 1;

    static constexpr double pressure(Mat<n_var> state)
    {
      double mmtm_sq = 0.;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        mmtm_sq += (state(j_dim))*(state(j_dim));
      }
      return (heat_rat - 1.)*((state(n_dim + 1)) - 0.5*mmtm_sq/(state(n_dim)));
    }

    static constexpr Mat<n_update> flux(Mat<n_var> state, Mat<n_dim> normal)
    {
      Mat<n_var> f;
      f(n_dim) = 0.;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        f(n_dim) += state(j_dim)*normal(j_dim);
      }
      double scaled = f(n_dim)/state(n_dim);
      double pres = pressure(state);
      HEXED_ASSERT_THERM_ADMIS
      f(n_var - 1) = (state(n_dim + 1) + pres)*scaled;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        f(j_dim) = state(j_dim)*scaled + pres*normal(j_dim);
      }
      return f;
    }

    static constexpr Mat<n_update> flux_num(Mat<n_var, 2> face_state, Mat<n_dim> normal)
    {
      Mat<n_var, 2> face_flux;
      Mat<2> vol_flux;
      Mat<2> sound_speed;
      double nrml_mag = normal.norm();
      for (int i_side = 0; i_side < 2; ++i_side) {
        auto f = face_flux(Eigen::all, i_side);
        auto state = face_state(Eigen::all, i_side);
        f(n_dim) = 0.;
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          f(n_dim) += state(j_dim)*normal(j_dim);
        }
        double scaled = f(n_dim)/state(n_dim);
        double pres = pressure(state);
        HEXED_ASSERT_THERM_ADMIS
        f(n_var - 1) = (state(n_dim + 1) + pres)*scaled;
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          f(j_dim) = state(j_dim)*scaled + pres*normal(j_dim);
        }
        vol_flux(i_side) = scaled;
        sound_speed(i_side) = std::sqrt(heat_rat*pres/state(n_dim))*nrml_mag;
      }
      Mat<2> wave_speed;
      wave_speed(0) = std::min(vol_flux(0) - sound_speed(0), vol_flux(1) - sound_speed(1));
      wave_speed(1) = std::max(vol_flux(0) + sound_speed(0), vol_flux(1) + sound_speed(1));
      return hll(wave_speed, face_flux, face_state);
    }

    static constexpr Mat<n_dim, n_update> flux_visc(Mat<n_dim, n_var> grad, Mat<n_dim, n_dim> nrmls, double av_coef)
    {
      return -av_coef*nrmls*grad;
    }
  };
};

// represents the nonuniform linear advection equation
// used in the smoothness-based artificial viscosity scheme
template <int n_dim>
class Advection
{
  public:
  Advection() = delete;
  static constexpr bool is_viscous = false;
  static constexpr bool has_convection = true;
  static constexpr int n_var = n_dim + 1;
  static constexpr int curr_start = n_dim;
  static constexpr int ref_start = n_dim + 1;
  static constexpr int n_update = 1;
  static constexpr int tss_pow = 1;

  static constexpr Mat<1> flux(Mat<n_var> state, Mat<n_dim> normal)
  {
    double nrml_veloc = 0;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      nrml_veloc += state(j_dim)*normal(j_dim);
    }
    return Mat<1>::Constant(nrml_veloc*state(n_dim));
  }

  static constexpr Mat<1> flux_num(Mat<n_var, 2> face_vars, Mat<n_dim> normal)
  {
    Mat<1, 2> vol_flux = normal.transpose()*face_vars(Eigen::seqN(0, n_dim), Eigen::all);
    Mat<1, 2> face_state = face_vars(n_dim, Eigen::all);
    Mat<1, 2> face_flux = face_state.cwiseProduct(vol_flux);
    Mat<2> wave_speed;
    for (int i_side = 0; i_side < 2; ++i_side) {
      int sign = 2*i_side - 1;
      wave_speed(i_side) = std::min(vol_flux(0) + sign*.01, vol_flux(1) + sign*.01);
    }
    return hll(wave_speed, face_flux, face_state);
  }
};

template <int n_dim>
class Fix_therm_admis
{
  public:
  Fix_therm_admis() = delete;
  static constexpr bool is_viscous = true;
  static constexpr bool has_convection = false;
  static constexpr int n_var = n_dim + 2;
  static constexpr int curr_start = 0;
  static constexpr int ref_start = n_var;
  static constexpr int visc_start = 2*n_var;
  static constexpr int n_update = n_var;
  static constexpr int tss_pow = 2;

  static constexpr Mat<n_update> flux_num(Mat<n_var, 2> face_state, Mat<n_dim> normal) {return Mat<n_update>::Zero();}
  static constexpr Mat<n_dim, n_update> flux_visc(Mat<n_dim, n_var> grad, Mat<n_dim, n_dim> nrmls, double av_coef)
  {
    return -nrmls*grad;
  }
};

}
#endif
