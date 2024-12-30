#ifndef HEXED_PDE_HPP_
#define HEXED_PDE_HPP_

#include "math.hpp"
#include "constants.hpp"
#include "Transport_model.hpp"
#include "Gauss_legendre.hpp"

/*! \brief This namespace contains classes representing the different PDEs Hexed can solve.
 * \details They are all possible arguments to the `Spatial` class template.
 * They define the organization of the state data, fluxes, and speeds of information
 * propagation for computing time steps.
 */
namespace hexed::pde {

constexpr int tss_offset(int n_var) {return n_var + 0;}
constexpr int bulk_av_offset(int n_var) {return n_var + 1;}
constexpr int laplacian_av_offset(int n_var) {return n_var + 2;}
constexpr int forcing_offset(int n_var) {return n_var + 3;}
constexpr int advection_offset(int n_var) {return n_var + 7;}

/*!
 * contains a PDE class representing the Naver-Stokes equations
 * with template options to specify the details of the equation set
 */
template <bool visc = false, Turbulence_model turb = laminar>
class Navier_stokes {
  public:
  Navier_stokes() = delete;

  template <int n_dim, int row_size>
  class Pde {
    const int _n_var;
    public:
    static constexpr bool has_diffusion = visc;
    static constexpr bool has_convection = true;
    static constexpr bool has_source = visc && (turb != laminar);
    static constexpr int n_update = n_dim + 2 + 2*(turb == k_omega);
    static constexpr int n_state = n_dim + 4 + 2*(turb == k_omega);
    static constexpr int n_extrap = n_dim + 2 + 2*(turb == k_omega);
    static constexpr int i_mass = n_dim;
    static constexpr int i_energy = n_dim + 1;
    static constexpr int i_turb_kin_ener = n_dim + 2;
    static constexpr int i_turb_diss = n_dim + 3;
    static constexpr int i_bulk_art_visc = n_update;
    static constexpr int i_laplacian_art_visc = n_update + 1;
    static constexpr double heat_rat = 1.4;

    static constexpr double alpha = 13./25.;
    static constexpr double beta_s = 9./100.;
    static constexpr double beta_0 = 0.0708;
    static constexpr double sigma = 1./2.;
    static constexpr double sigma_s = 3./5.;
    static constexpr double sigma_do = 1./8.;
    static constexpr double c_lim = 7./8.;
    static constexpr double turb_prandtl = .9;

    Transport_model dyn_visc;
    Transport_model therm_cond;

    Pde(int n_var, Transport_model dynamic_visc = inviscid, Transport_model thermal_cond = inviscid)
    : _n_var{n_var}, dyn_visc{dynamic_visc}, therm_cond{thermal_cond}
    {}

    Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
      Mat<n_extrap> extrap;
      for (int i_var = 0; i_var < n_extrap; ++i_var) extrap(i_var) = data[i_var*stride];
      return extrap;
    }

    void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
      for (int i_var = 0; i_var < n_update; ++i_var) data[i_var*stride] += update(i_var);
    }

    /*! \todo __Carter:__ This is the class you have to modify to implement the \f$ k\text{-}\omega \f$ equations.
     * Currently, to test that the numerical kernel is capable of handling the extra equations,
     * I've set it to solve the arbitrary equations:
     * \f$
     * \frac{\partial \rho k}{\partial t}
     * + \frac{\partial}{\partial x_j}
     *   \left( \rho u_j k - \frac{\mu}{\rho} \frac{\partial \rho k}{\partial x_j} \right)
     * = -0.1 \mu k
     * \f$
     * and
     * \f$
     * \frac{\partial \rho \tilde{\omega}}{\partial t}
     * + \frac{\partial}{\partial x_j}
     *   \left( \rho u_j \tilde{\omega} - \frac{\mu}{\rho} \frac{\partial \rho \tilde{\omega}}{\partial x_j} \right)
     * = -0.1 \mu \tilde{\omega}
     * \f$
     * (Einstein summation convention)
     */
    template <int n_dim_flux>
    class Computation {
      const Pde& _eq;
      public:
      Mat<config::debug_variables> debug_variables; //!< \brief can be populated at any time
      Computation(const Pde& eq) : _eq{eq} {}

      Mat<n_state> state;
      void fetch_state(int stride, const double* data) {
        for (int i_var = 0; i_var < n_state - 2; ++i_var) state(i_var) = data[i_var*stride];
        state(i_bulk_art_visc) = data[bulk_av_offset(_eq._n_var)*stride];
        state(i_laplacian_art_visc) = data[laplacian_av_offset(_eq._n_var)*stride];
      }
      Mat<n_update> update_state;
      void fetch_extrap_state(int stride, const double* data) {
        state.setZero();
        for (int i_var = 0; i_var < n_extrap; ++i_var) state(i_var) = data[i_var*stride];
        update_state = state(Eigen::seqN(0, n_update));
      }

      double mass;
      double kin_ener;
      double pressure;
      void compute_scalars_conv() {
        mass = state(i_mass);
        kin_ener = 0;
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) kin_ener += state(i_dim)*state(i_dim);
        kin_ener *= .5/mass;
        pressure = (heat_rat - 1.)*(state(i_energy) - kin_ener);
      }

      Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
      Mat<n_update, n_dim_flux> flux_conv;
      void compute_flux_conv() {
        compute_scalars_conv();
        for (int i_dim = 0; i_dim < n_dim_flux; ++i_dim) {
          flux_conv(i_mass, i_dim) = 0;
          for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
            flux_conv(i_mass, i_dim) += state(j_dim)*normal(j_dim, i_dim);
          }
          double vol_flux = flux_conv(i_mass, i_dim)/mass;
          flux_conv(i_energy, i_dim) = (state(i_energy) + pressure)*vol_flux;
          if constexpr (turb == k_omega) {
            // these convective fluxes for the turbulent variables should be correct
            flux_conv(i_turb_kin_ener, i_dim) = state(i_turb_kin_ener)*vol_flux;
            flux_conv(i_turb_diss, i_dim) = state(i_turb_diss)*vol_flux;
          }
          for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
            flux_conv(j_dim, i_dim) = state(j_dim)*vol_flux + pressure*normal(j_dim, i_dim);
          }
        }
      }

      double bulk_av;
      double laplacian_av;
      double sqrt_temp;
      double dyn_visc_coef;
      double therm_cond_coef;
      double energy_cond;
      double real_turb_diss;
      double k_bar;
      double beta;
      //! \todo __Carter:__ Compute whatever variables you need for the turbulent fluxes
      //! which might also be needed for source terms and/or the time step calculation.
      void compute_scalars_diff() {
        bulk_av = std::abs(state(i_bulk_art_visc));
        laplacian_av = std::abs(state(i_laplacian_art_visc));
        double spec_heat_v = constants::specific_gas_air/(heat_rat - 1.);
        double spec_heat_p = heat_rat*constants::specific_gas_air;
        sqrt_temp = std::sqrt(std::max((state(i_energy) - kin_ener)/mass, 0.)/spec_heat_v);
        // taking abs ensures that they will never be negative
        // and makes the probability that they are exactly 0 very low, which is good cause we have to divide by them
        real_turb_diss = std::max(10., std::exp(state(i_turb_diss)/mass));
        k_bar = std::abs(state(i_turb_kin_ener)/mass);
        dyn_visc_coef = _eq.dyn_visc.coefficient(sqrt_temp);
        therm_cond_coef = _eq.therm_cond.coefficient(sqrt_temp);
        energy_cond = therm_cond_coef/spec_heat_v;
      }

      Mat<n_extrap, n_dim> gradient;
      Mat<n_update, n_dim_flux> flux_diff;
      /*! \todo __Carter:__ modify `flux_diff` to include turbulence modeling.
       * Set `flux_diff_phys(i_turb_kin_ener)` and `flux_diff_phys(i_turb_diss)` to contain the source terms of
       * \f$ \rho k \f$ and \f$ \rho \tilde{\omega} \f$, respectively.
       * Also modify the other fluxes to include turbulent terms.
       */
      void compute_flux_diff() {
        compute_scalars_diff();
        auto seq = Eigen::seqN(0, n_dim);
        auto mmtm = state(seq);
        Mat<n_dim> veloc = mmtm/mass;
        Mat<n_dim, n_dim> veloc_grad = (gradient(seq, all) - veloc*gradient(i_mass, all))/mass;
        double divergence = veloc_grad.trace();
        Mat<n_dim, n_dim> rotation = .5*(veloc_grad - veloc_grad.transpose());
        Mat<n_dim, n_dim> strain_rate = .5*(veloc_grad + veloc_grad.transpose());
        Mat<n_dim, n_dim> identity = Mat<n_dim, n_dim>::Identity();

        /*! \note instead of \f$ \hat{\omega} = \max(\omega, C_{lim} \sqrt{2 \bar{S}_{ij} \bar{S}_{ij}/\beta^*}) \f$,
         * do \f$ \hat{\omega} = (\omega^4 + (C_{lim} \sqrt{2 \bar{S}_{ij} \bar{S}_{ij}/\beta^*})^4)^{1/4} \f$
         * for a smoother transition.
         */
        double lim_sq = c_lim*c_lim*2*(strain_rate - 1./3.*divergence*identity).squaredNorm()/beta_s;
        double omega_hat = std::pow(math::pow(real_turb_diss, 4) + lim_sq*lim_sq, .25);

        Mat<n_dim, n_dim> turb_stress_per_k = 2*mass/omega_hat*strain_rate - 2./3.*mass*(1. + divergence/omega_hat)*identity;
        Mat<n_dim, n_dim> stress = 2*dyn_visc_coef*strain_rate + (bulk_av*mass - 2./3.*dyn_visc_coef)*divergence*identity
                                   + turb_stress_per_k*k_bar;
        Mat<n_update, n_dim> flux_diff_phys = -laplacian_av*gradient; // flux in physical space
        flux_diff_phys(seq, all) -= stress;
        Mat<1, n_dim> int_ener_grad = -state(i_energy)/mass/mass*gradient(i_mass, all)
                                      + gradient(i_energy, all)/mass - veloc.transpose()*veloc_grad;
        flux_diff_phys(i_energy, all) -= veloc.transpose()*stress
                                         + (energy_cond + heat_rat*mass*k_bar/omega_hat/turb_prandtl)*int_ener_grad;
        source.setZero();
        if constexpr (turb == k_omega) {
          Mat<n_dim> grad_k = gradient(i_turb_kin_ener, all)/mass
                              - state(i_turb_kin_ener)/mass/mass*gradient(i_mass, all);
          Mat<n_dim> grad_omega = gradient(i_turb_diss, all)/mass
                                  - state(i_turb_diss)/mass/mass*gradient(i_mass, all);
          // note that unlike in the turbulent viscosity, here we do _not_ use `omega_hat`
          flux_diff_phys(i_turb_kin_ener, all) -= (dyn_visc_coef + sigma_s*mass*k_bar/real_turb_diss)*grad_k;
          flux_diff_phys(i_turb_diss, all) -= (dyn_visc_coef + sigma*mass*k_bar/real_turb_diss)*grad_omega;

          Mat<n_dim, n_dim> S_hat = strain_rate - .5*divergence*identity;
          double sum = 0;
          for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
              for (int k = 0; k < n_dim; ++k){
                sum += rotation(i, j)*rotation(j, k)*S_hat(k, i);
              }
            }
          }
          double chi_o = std::abs(sum)/math::pow(beta_s*real_turb_diss, 3);
          double f_beta = (1. + 85.*chi_o)/(1. + 100.*chi_o);
          beta = beta_0*f_beta;

          double prod_per_k = 0.0;
          for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
              prod_per_k += turb_stress_per_k(i, j)*veloc_grad(i, j);
            }
          }
          double lim = 1e4*mass*real_turb_diss;
          double lim_factor = lim/std::sqrt(lim*lim + prod_per_k*prod_per_k);
          debug_variables(0) = mass*k_bar/omega_hat;
          prod_per_k = lim_factor*prod_per_k;
          double dot = grad_k.dot(grad_omega);
          double grad_k_omega_source = std::max(sigma_do*mass/real_turb_diss*grad_k.dot(grad_omega), 0.);
          double grad_omega_source = (dyn_visc_coef + sigma*mass*k_bar/real_turb_diss)*grad_omega.squaredNorm();
          source(i_turb_kin_ener) = prod_per_k*k_bar - beta_s*real_turb_diss*state(i_turb_kin_ener);
          source(i_turb_diss) = alpha*prod_per_k + grad_omega_source + grad_k_omega_source - beta*mass*real_turb_diss;
          source(i_energy) = -source(i_turb_kin_ener);
        }
        flux_diff = flux_diff_phys*normal; // flux in reference space
      }

      Mat<n_update> source;
      //! \brief does nothing---source terms are computed in `compute_diffusion`
      //! \details because they need access to gradients
      void compute_source() {}

      double char_speed;
      void compute_char_speed() {
        // numerical estimate (not less than actual speed of sound)
        const double sound_speed = std::sqrt(heat_rat*(heat_rat - 1)*state(i_energy)/state(i_mass));
        const double speed = state(Eigen::seqN(0, n_dim)).norm()/state(i_mass);
        char_speed = sound_speed + speed;
      }

      double diffusivity;
      void compute_diffusivity() {
        compute_scalars_conv();
        compute_scalars_diff();
        // this is a conservative estimate for `dyn_visc_turb` because `omega_hat` is not available
        double dyn_visc_turb = (turb == k_omega) ? std::abs(state(i_turb_kin_ener))*std::exp(-state(i_turb_diss)/state(i_mass)) : 0.;
        diffusivity = 10*std::abs(laplacian_av) + math::max(
          (dyn_visc_coef + math::max(1, sigma, sigma_s)*dyn_visc_turb)/mass,
          std::abs(bulk_av) + (dyn_visc_coef + dyn_visc_turb)/mass,
          (dyn_visc_coef + dyn_visc_turb)/mass + (energy_cond + heat_rat*dyn_visc_turb/turb_prandtl)/mass
        );
      }

      double decay;
      void compute_decay() {
        decay = 0;
        if constexpr (turb == k_omega) {
          real_turb_diss = std::exp(state(i_turb_diss)/mass);
          decay = std::max(beta_s*real_turb_diss, beta_s); // note: beta <= beta_s
        }
      }
    };

    /*!
     * Decomposes state vectors into characteristics
     * which are eigenvectors of the Jacobian of the inviscid flux function.
     * This is useful for characteristic-based boundary conditions.
     */
    class Characteristics {
      static constexpr int n_var_euler = n_dim + 2;
      Mat<3> vals;
      Mat<3, 3> vecs;
      // using a QR factorization allows a least-squares solution to be found if matrix is singular (i.e. if pressure is 0)
      Eigen::ColPivHouseholderQR<Mat<3, 3>> fact;
      Mat<n_dim> dir; // normalized flux direction
      // some properties of the reference state
      double mass;
      Mat<n_dim> veloc;
      double nrml(Mat<n_dim> vec) {return dir.dot(vec);}
      Mat<n_dim> tang(Mat<n_dim> vec) {return vec - dir*nrml(vec);}

      public:
      /*!
       * construct with a direction in which to compute the flux
       * and a reference state vector about which to compute the Jacobian
       */
      Characteristics(Mat<n_var_euler> state, Mat<n_dim> direction)
      : dir{direction/direction.norm()}
      , mass{state(i_mass)}
      , veloc{state(Eigen::seqN(0, n_dim))/mass}
      {
        // compute more properties of the reference state
        double vsq = veloc.squaredNorm();
        double pres = .4*(state(i_energy) - .5*mass*vsq);
        double sound_speed = std::sqrt(1.4*std::max(pres, 0.)/mass);
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
        fact.compute(vecs);
        HEXED_ASSERT(fact.info() == Eigen::Success, "QR factorization failed");
      }
      //! \brief get eigenvalues of Jacobian
      inline Mat<3> eigvals() {return vals;}
      /*! \brief Decompose a state vector into eigenspaces.
       * Column `j` should be an eigenvector of the Jacobian with eigenvalue `eigvals()(j)`
       * and the sum of the columns should be `state`.
       */
      Mat<n_var_euler, 3> decomp(Mat<n_var_euler> state) {
        Mat<n_dim> mmtm = state(Eigen::seqN(0, n_dim));
        // component of tangential momentum perturbation which is not induced by mass perturbation
        Mat<n_dim> mmtm_correction = tang(mmtm) - state(i_mass)*tang(veloc);
        // compute state for 1D eigenvector problem
        Mat<3> state_1d;
        state_1d <<
          nrml(mmtm),
          state(i_mass),
          state(i_energy) - veloc.dot(mmtm_correction); //! \todo shouldn't this have a `0.5*`?
        // decompose 1D state into eigenvectors
        Mat<1, 3> eig_basis = fact.solve(state_1d).transpose();
        Mat<3, 3> eig_decomp = vecs.array().rowwise()*eig_basis.array();
        // ND eigenvector decomposition
        Mat<n_var_euler, 3> d(state.rows(), 3);
        d(Eigen::seqN(n_dim, 2), Eigen::all) = eig_decomp(Eigen::seqN(1, 2), Eigen::all);
        d(Eigen::seqN(0, n_dim), Eigen::all) = dir*eig_decomp(0, Eigen::all) + tang(veloc)*eig_basis; // second term accounts for tangential momentum induced by mass perturbation
        d(Eigen::seqN(0, n_dim), 2) += mmtm_correction;
        d(i_energy, 2) += veloc.dot(mmtm_correction); // correct energy to account for tangential momentum perturbation
        return d;
      }
    };
  };
};

/*!
 * represents the nonuniform linear advection equation
 * used in the smoothness-based artificial viscosity scheme
 */
template <int n_dim, int row_size>
class Advection {
  const int _n_var;
  static constexpr int _n_adv = row_size;
  const double _advect_length;
  Mat<row_size> _nodes;

  public:
  static constexpr bool has_diffusion = false;
  static constexpr bool has_convection = true;
  static constexpr bool has_source = true;
  static constexpr int n_state = n_dim + _n_adv;
  static constexpr int n_extrap = n_dim + _n_adv;
  static constexpr int n_update = _n_adv;

  Advection(int n_var, double advect_length)
  : _n_var{n_var}, _advect_length{advect_length}, _nodes{2*Gauss_legendre(row_size).nodes() - Mat<row_size>::Ones()}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_dim; ++i_var) extrap(i_var) = data[i_var*stride];
    for (int i_adv = 0; i_adv < _n_adv; ++i_adv) extrap(n_dim + i_adv) = data[(advection_offset(_n_var) + i_adv)*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool is_critical) const {
    double pseudo = 1 + data[tss_offset(_n_var)*stride]*2/_advect_length;
    for (int i_adv = 0; i_adv < _n_adv; ++i_adv) {
      double& d = data[(advection_offset(_n_var) + i_adv)*stride];
      if (is_critical) d = (d + update(i_adv))/pseudo;
      else d += update(i_adv)/pseudo;
    }
  }

  template <int n_dim_flux>
  class Computation {
    const Advection& _eq;
    public:
    Mat<config::debug_variables> debug_variables; //!< \brief can be populated at any time
    Computation(const Advection& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {
      state = _eq.fetch_extrap(stride, data);
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_extrap; ++i_var) state(i_var) = data[i_var*stride];
      for (int i_adv = 0; i_adv < _n_adv; ++i_adv) update_state(i_adv) = state(n_dim + i_adv);
    }

    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_update, n_dim_flux> flux_conv;
    void compute_flux_conv() {
      for (int i_dim = 0; i_dim < n_dim_flux; ++i_dim) {
        double nrml_veloc = 0;
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          nrml_veloc += state(j_dim)*normal(j_dim, i_dim);
        }
        for (int i_adv = 0; i_adv < _n_adv; ++i_adv) {
          flux_conv(i_adv, i_dim) = _eq._nodes(i_adv)*nrml_veloc*state(n_dim + i_adv);
        }
      }
    }

    Mat<n_update> source;
    void compute_source() {
      source.setConstant(2/_eq._advect_length);
    }

    double char_speed;
    void compute_char_speed() {
      char_speed = std::max(1., state(Eigen::seqN(0, n_dim)).norm());
    }

    double decay;
    void compute_decay() {
      decay = 0;
    }
  };
};

/*!
 * represents the uniform linear diffusion equation
 * used in the smoothness-based artificial viscosity scheme
 */
template <int n_dim, int row_size>
class Smooth_art_visc {
  const int _n_var;
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = false;
  static constexpr bool has_source = true;
  static constexpr int n_state = 4;
  static constexpr int n_extrap = 3;
  static constexpr int n_update = 3;
  const double _diff_time;
  const double _cheby;

  Smooth_art_visc(int n_var, double diff_time, double chebyshev_step)
  : _n_var{n_var}, _diff_time{diff_time}, _cheby{chebyshev_step}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_extrap; ++i_var) extrap(i_var) = data[(forcing_offset(_n_var) + 1 + i_var)*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
    double pseudo = 1 + data[tss_offset(_n_var)*stride]*_cheby/_diff_time;
    for (int i_var = 0; i_var < n_update; ++i_var) {
      double& d = data[(forcing_offset(_n_var) + 1 + i_var)*stride];
      d += update(i_var);
      if (critical) d /= pseudo;
    }
  }

  template <int n_dim_flux>
  class Computation {
    const Smooth_art_visc& _eq;
    public:
    Mat<config::debug_variables> debug_variables; //!< \brief can be populated at any time
    Computation(const Smooth_art_visc& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_state; ++i_var) state(i_var) = data[(forcing_offset(_eq._n_var) + i_var)*stride];
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_update; ++i_var) update_state(i_var) = data[i_var*stride];
    }

    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_extrap, n_dim> gradient;
    Mat<n_update, n_dim_flux> flux_diff;
    void compute_flux_diff() {
      flux_diff.noalias() = -gradient*normal;
    }

    double diffusivity;
    void compute_diffusivity() {
      diffusivity = 1;
    }

    Mat<n_update> source;
    void compute_source() {
      for (int i_var = 0; i_var < n_update; ++i_var) {
        double f = std::abs(state(i_var));
        source(i_var) = ((i_var == 1) ? std::sqrt(f) : f)/_eq._diff_time;
      }
    }

    double decay;
    void compute_decay() {
      decay = 0;
    }
  };
};

/*!
 * represents the uniform linear diffusion equation
 * used for fixing thermodynamic admissibility
 */
template <int n_dim, int row_size>
class Fix_therm_admis {
  int _n_var;
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = false;
  static constexpr bool has_source = false;
  static constexpr int n_state = n_dim + 2;
  static constexpr int n_update = n_state;
  static constexpr int n_extrap = n_state;

  Fix_therm_admis(int n_var) : _n_var{n_var} {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_extrap; ++i_var) extrap(i_var) = data[i_var*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
    for (int i_var = 0; i_var < n_update; ++i_var) data[i_var*stride] += update(i_var);
  }

  template <int n_dim_flux>
  class Computation {
    const Fix_therm_admis& _eq;
    public:
    Mat<config::debug_variables> debug_variables; //!< \brief can be populated at any time
    Computation(const Fix_therm_admis& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_state; ++i_var) state(i_var) = data[i_var*stride];
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_state; ++i_var) update_state(i_var) = data[i_var*stride];
    }

    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_extrap, n_dim> gradient;
    Mat<n_update, n_dim_flux> flux_diff;
    void compute_flux_diff() {
      flux_diff.noalias() = -gradient*normal;
    }

    double diffusivity;
    void compute_diffusivity() {
      diffusivity = 1;
    }
  };
};

}
#endif
