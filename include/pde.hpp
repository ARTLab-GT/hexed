#ifndef HEXED_PDE_HPP_
#define HEXED_PDE_HPP_

#include <memory>
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

class Transport_model
{
  public:
  virtual double operator()(double temperature) const = 0;
  virtual ~Transport_model() = default;
};

class Const_transport : public Transport_model
{
  double v;
  public:
  inline Const_transport(double value) : v{value} {}
  virtual inline double operator()(double temperature) const {return v;}
};

class Sutherland : public Transport_model
{
  double ref_val;
  double ref_temp;
  double offset;
  public:
  inline Sutherland(double reference_value, double reference_temperature, double temperature_offset)
  : ref_val{reference_value}, ref_temp{reference_temperature}, offset{temperature_offset}
  {}
  virtual inline double operator()(double temperature) const
  {
    return ref_val*std::pow(temperature/ref_temp, 1.5)*(ref_temp + offset)/(temperature + offset);
  }
};

std::shared_ptr<Transport_model> inviscid();
std::shared_ptr<Transport_model> air_const_dyn_visc();
std::shared_ptr<Transport_model> air_const_therm_cond();
std::shared_ptr<Transport_model> air_sutherland_dyn_visc();
std::shared_ptr<Transport_model> air_sutherland_therm_cond();

// contains a PDE class representing the Naver-Stokes equations
// with template options to specify the details of the equation set
template <bool visc = false>
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
    static constexpr bool is_viscous = visc;
    static constexpr bool has_convection = true;
    static constexpr int n_var = n_dim + 2;
    static constexpr int curr_start = 0;
    static constexpr int ref_start = n_var;
    static constexpr int visc_start = 2*n_var;
    static constexpr int n_update = n_var;
    static constexpr double heat_rat = 1.4;
    static constexpr double gas_const = 287.0528;
    static constexpr int tss_pow = 1;
    std::shared_ptr<Transport_model> dyn_visc;
    std::shared_ptr<Transport_model> therm_cond;

    Pde(std::shared_ptr<Transport_model> dyn_visc_model = inviscid(), std::shared_ptr<Transport_model> therm_cond_model = inviscid())
    : dyn_visc{dyn_visc_model},
      therm_cond{therm_cond_model}
    {}

    double pressure(Mat<n_var> state) const
    {
      double mmtm_sq = 0.;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        mmtm_sq += (state(j_dim))*(state(j_dim));
      }
      return (heat_rat - 1.)*((state(n_dim + 1)) - 0.5*mmtm_sq/(state(n_dim)));
    }

    // compute the convective flux
    Mat<n_update> flux(Mat<n_var> state, Mat<n_dim> normal) const
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
    Mat<n_update> flux_num(Mat<n_var, 2> face_state, Mat<n_dim> normal) const
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
    Mat<n_dim, n_update> flux_visc(Mat<n_var> state, Mat<n_dim, n_var> grad, double av_coef) const
    {
      auto seq = Eigen::seqN(0, n_dim);
      auto all = Eigen::all;
      auto mmtm = state(seq);
      double mass = state(n_dim);
      Mat<n_dim> veloc = mmtm/mass;
      Mat<n_dim, n_dim> veloc_grad = (grad(all, seq) - grad(all, n_dim)*veloc.transpose())/mass;
      double temp = (state(n_dim + 1)/mass - .5*veloc.squaredNorm())*(heat_rat - 1)/gas_const;
      Mat<n_dim, n_dim> stress = (*dyn_visc)(temp)*(veloc_grad + veloc_grad.transpose() - 2./3.*veloc_grad.trace()*Mat<n_dim, n_dim>::Identity());
      Mat<n_dim, n_update> flux = -av_coef*grad;
      flux(all, seq) -= stress;
      Mat<n_dim> temp_grad = -(state(n_dim + 1)/mass/mass*grad(all, n_dim) + grad(all, n_dim + 1)/mass - veloc_grad*veloc)/(gas_const/(heat_rat - 1));
      flux(all, n_dim + 1) -= stress*veloc + (*therm_cond)(temp)*temp_grad;
      return flux;
    }

    // maximum characteristic speed for convection
    double char_speed(Mat<n_var> state) const
    {
      const double sound_speed = std::sqrt(heat_rat*pressure(state)/state(n_dim));
      auto mmtm = state(Eigen::seqN(0, n_dim));
      const double veloc = std::sqrt(mmtm.dot(mmtm)/state(n_dim)/state(n_dim));
      return sound_speed + veloc;
    }

    // maximum effective diffusivity (for enforcing the CFL condition)
    double diffusivity(Mat<n_var> state, double av_coef) const
    {
      double mass = state(n_dim);
      Mat<n_dim> veloc = state(Eigen::seqN(0, n_dim))/mass;
      double temp = (state(n_dim + 1)/mass - .5*veloc.squaredNorm())*(heat_rat - 1)/gas_const;
      return av_coef + std::max((*dyn_visc)(temp)/mass, (*therm_cond)(temp)/(mass*gas_const/(heat_rat - 1)));
    }

    /*
     * Decomposes state vectors into characteristics
     * which are eigenvectors of the Jacobian of the inviscid flux function.
     * This is useful for characteristic-based boundary conditions.
     */
    class Characteristics
    {
      Mat<3> vals;
      Mat<3, 3> vecs;
      Mat<3, 3> vecs_inv;
      Mat<n_dim> dir; // normalized flux direction
      // some properties of the reference state
      double mass;
      Mat<n_dim> veloc;
      double nrml(Mat<n_dim> vec) {return dir.dot(vec);}
      Mat<n_dim> tang(Mat<n_dim> vec) {return vec - dir*nrml(vec);}

      public:
      // construct with a direction in which to compute the flux
      // and a reference state vector about which to compute the Jacobian
      Characteristics(Mat<n_var> state, Mat<n_dim> direction)
      : dir{direction/direction.norm()},
        mass{state(n_dim)},
        veloc{state(Eigen::seqN(0, n_dim))/mass}
      {
        // compute more properties of the reference state
        double vsq = veloc.squaredNorm();
        double pres = .4*(state(n_dim + 1) - .5*mass*vsq);
        double sound_speed = std::sqrt(1.4*pres/mass);
        // compute eigenvalues
        vals(2) = nrml(veloc);
        vals(0) = vals(2) - sound_speed;
        vals(1) = vals(2) + sound_speed;
        // compute 1D eigenvectors
        double d_mass = 1;
        for (int sign = 0; sign < 2; ++sign) {
          double d_veloc = (2*sign - 1)*sound_speed/mass*d_mass;
          double d_pres = 1.4*pres/mass*d_mass;
          vecs(Eigen::all, sign) <<
            d_mass*vals(2) + mass*d_veloc,
            d_mass,
            d_pres/.4 + .5*d_mass*vsq + mass*vals(2)*d_veloc;
        }
        vecs(Eigen::all, 2) << d_mass*vals(2), d_mass, .5*d_mass*vsq;
        vecs_inv = vecs.inverse();
      }
      // get eigenvalues of Jacobian
      inline Mat<3> eigvals() {return vals;}
      /*
       * Decompose a state vector into eigenspaces.
       * Column `j` should be an eigenvector of the Jacobian with eigenvalue `eigvals()(j)`
       * and the sum of the columns should be `state`.
       */
      Mat<n_var, 3> decomp(Mat<n_var> state)
      {
        Mat<n_dim> mmtm = state(Eigen::seqN(0, n_dim));
        // component of tangential momentum perturbation which is not induced by mass perturbation
        Mat<n_dim> mmtm_correction = tang(mmtm) - state(n_dim)*tang(veloc);
        // compute state for 1D eigenvector problem
        Mat<3> state_1d;
        state_1d <<
          nrml(mmtm),
          state(n_dim),
          state(n_dim + 1) - veloc.dot(mmtm_correction);
        // decompose 1D state into eigenvectors
        Mat<1, 3> eig_basis = (vecs_inv*state_1d).transpose();
        Mat<3, 3> eig_decomp = vecs.array().rowwise()*eig_basis.array();
        // ND eigenvector decomposition
        Mat<n_var, 3> d(state.rows(), 3);
        d(Eigen::seqN(n_dim, 2), Eigen::all) = eig_decomp(Eigen::seqN(1, 2), Eigen::all);
        d(Eigen::seqN(0, n_dim), Eigen::all) = dir*eig_decomp(0, Eigen::all) + tang(veloc)*eig_basis; // second term accounts for tangential momentum induced by mass perturbation
        d(Eigen::seqN(0, n_dim), 2) += mmtm_correction;
        d(n_dim + 1, 2) += veloc.dot(mmtm_correction); // correct energy to account for tangential momentum perturbation
        return d;
      }
    };
  };
  #undef ASSERT_THERM_ADMIS
};

// represents the nonuniform linear advection equation
// used in the smoothness-based artificial viscosity scheme
template <int n_dim>
class Advection
{
  public:
  static constexpr bool is_viscous = false;
  static constexpr bool has_convection = true;
  static constexpr int n_var = n_dim + 1;
  static constexpr int curr_start = n_dim;
  static constexpr int ref_start = n_dim + 1;
  static constexpr int n_update = 1;
  static constexpr int tss_pow = 1;

  Mat<1> flux(Mat<n_var> state, Mat<n_dim> normal) const
  {
    double nrml_veloc = 0;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      nrml_veloc += state(j_dim)*normal(j_dim);
    }
    return Mat<1>::Constant(nrml_veloc*state(n_dim));
  }

  Mat<1> flux_num(Mat<n_var, 2> face_vars, Mat<n_dim> normal) const
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

  double char_speed(Mat<n_var> state) const
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
  static constexpr bool is_viscous = true;
  static constexpr bool has_convection = false;
  static constexpr int n_var = 1;
  static constexpr int curr_start = 0;
  static constexpr int ref_start = 1;
  static constexpr int visc_start = 2;
  static constexpr int n_update = 1;
  static constexpr int tss_pow = 2;

  Mat<n_update> flux_num(Mat<n_var, 2> face_state, Mat<n_dim> normal) const {return Mat<n_update>::Zero();}
  Mat<n_dim, n_update> flux_visc(Mat<n_var> state, Mat<n_dim, n_var> grad, double av_coef) const
  {
    return -grad;
  }
  double diffusivity(Mat<n_var> state, double av_coef) const {return 1;}
};

// represents the uniform linear diffusion equation
// used for fixing thermodynamic admissibility
template <int n_dim>
class Fix_therm_admis
{
  public:
  static constexpr bool is_viscous = true;
  static constexpr bool has_convection = false;
  static constexpr int n_var = n_dim + 2;
  static constexpr int curr_start = 0;
  static constexpr int ref_start = n_var;
  static constexpr int visc_start = 2*n_var;
  static constexpr int n_update = n_var;
  static constexpr int tss_pow = 2;

  Mat<n_update> flux_num(Mat<n_var, 2> face_state, Mat<n_dim> normal) const {return Mat<n_update>::Zero();}
  Mat<n_dim, n_update> flux_visc(Mat<n_var> state, Mat<n_dim, n_var> grad, double av_coef) const
  {
    return -grad;
  }
  double diffusivity(Mat<n_var> state, double av_coef) const {return 1;}
};

}
#endif
