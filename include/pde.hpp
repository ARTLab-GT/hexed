#ifndef HEXED_PDE_HPP_
#define HEXED_PDE_HPP_

#include "math.hpp"

/*
 * This namespace contains classes representing the different PDEs Hexed can solve.
 * They are all possible arguments to the `Spatial` class template.
 * They define the organization of the state data, fluxes, and speeds of information
 * propagation for computing time steps.
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
  // check that the flow state is thermodynamically admissible
  #define ASSERT_THERM_ADMIS \
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

    // compute the convective flux
    static constexpr Mat<n_update> flux(Mat<n_var> state, Mat<n_dim> normal)
    {
      Mat<n_var> f;
      f(n_dim) = 0.;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        f(n_dim) += state(j_dim)*normal(j_dim);
      }
      double scaled = f(n_dim)/state(n_dim);
      double pres = pressure(state);
      ASSERT_THERM_ADMIS
      f(n_var - 1) = (state(n_dim + 1) + pres)*scaled;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        f(j_dim) = state(j_dim)*scaled + pres*normal(j_dim);
      }
      return f;
    }

    // compute the numerical convective flux shared at the element faces
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
        ASSERT_THERM_ADMIS
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

    // compute the viscous flux
    static constexpr Mat<n_dim, n_update> flux_visc(Mat<n_var> state, Mat<n_dim, n_var> grad, double av_coef)
    {
      return -av_coef*grad;
    }

    // maximum characteristic speed for convection
    static constexpr double char_speed(Mat<n_var> state)
    {
      const double sound_speed = std::sqrt(heat_rat*pressure(state)/state(n_dim));
      auto mmtm = state(Eigen::seqN(0, n_dim));
      const double veloc = std::sqrt(mmtm.dot(mmtm)/state(n_dim)/state(n_dim));
      return sound_speed + veloc;
    }

    // maximum effective diffusivity (for enforcing the CFL condition)
    static constexpr double diffusivity(Mat<n_var> state, double av_coef)
    {
      return av_coef;
    }
  };
  #undef ASSERT_THERM_ADMIS
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

  static constexpr double char_speed(Mat<n_var> state)
  {
    return state(Eigen::seqN(0, n_dim)).norm();
  }
};

// represents the uniform linear diffusion equation
// used in the smoothness-based artificial viscosity scheme
template <int n_dim>
class Smooth_art_visc
{
  public:
  Smooth_art_visc() = delete;
  static constexpr bool is_viscous = true;
  static constexpr bool has_convection = false;
  static constexpr int n_var = 1;
  static constexpr int curr_start = 0;
  static constexpr int ref_start = 1;
  static constexpr int visc_start = 2;
  static constexpr int n_update = 1;
  static constexpr int tss_pow = 2;

  static constexpr Mat<n_update> flux_num(Mat<n_var, 2> face_state, Mat<n_dim> normal) {return Mat<n_update>::Zero();}
  static constexpr Mat<n_dim, n_update> flux_visc(Mat<n_var> state, Mat<n_dim, n_var> grad, double av_coef)
  {
    return -grad;
  }
  static constexpr double diffusivity(Mat<n_var> state, double av_coef) {return 1;}
};

// represents the uniform linear diffusion equation
// used for fixing thermodynamic admissibility
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
  static constexpr Mat<n_dim, n_update> flux_visc(Mat<n_var> state, Mat<n_dim, n_var> grad, double av_coef)
  {
    return -grad;
  }
  static constexpr double diffusivity(Mat<n_var> state, double av_coef) {return 1;}
};

}
#endif
