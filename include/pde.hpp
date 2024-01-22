#ifndef HEXED_PDE_HPP_
#define HEXED_PDE_HPP_

#include "math.hpp"
#include "constants.hpp"
#include "Transport_model.hpp"

/*! \brief This namespace contains classes representing the different PDEs Hexed can solve.
 * \details They are all possible arguments to the `Spatial` class template.
 * They define the organization of the state data, fluxes, and speeds of information
 * propagation for computing time steps.
 */
namespace hexed::pde
{

/*!
 * contains a PDE class representing the Naver-Stokes equations
 * with template options to specify the details of the equation set
 */
template <bool visc = false>
class Navier_stokes
{
  //! check that the flow state is thermodynamically admissible
  #define ASSERT_THERM_ADMIS \
    HEXED_ASSERT(state(n_dim) > 0, "nonpositive density", assert::Numerical_exception); \
    HEXED_ASSERT(state(n_dim + 1) >= 0, "negative energy", assert::Numerical_exception); \
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) { \
      HEXED_ASSERT(!std::isnan(state(n_dim)), "momentum is NaN", assert::Numerical_exception); \
    } \

  public:
  Navier_stokes() = delete;

  template <int n_dim>
  class Pde
  {
    bool _laplacian;
    public:
    static constexpr bool has_diffusion = visc;
    static constexpr bool has_convection = true;
    static constexpr bool has_source = false;
    static constexpr int n_update = n_dim + 2;
    static constexpr int n_state = n_dim + 3;
    static constexpr int n_extrap = n_dim + 2;
    static constexpr double heat_rat = 1.4;
    Transport_model dyn_visc;
    Transport_model therm_cond;

    Pde(Transport_model dynamic_visc = inviscid, Transport_model thermal_cond = inviscid, bool laplacian = false)
    : _laplacian{laplacian}, dyn_visc{dynamic_visc}, therm_cond{thermal_cond}
    {}

    Mat<n_extrap> fetch_extrap(int stride, const double* data) const
    {
      Mat<n_extrap> extrap;
      for (int i_var = 0; i_var < n_extrap; ++i_var) extrap(i_var) = data[i_var*stride];
      return extrap;
    }

    void write_update(Mat<n_update> update, int stride, double* data, bool critical) const
    {
      for (int i_var = 0; i_var < n_update; ++i_var) data[i_var*stride] += update(i_var);
    }

    template <int n_dim_flux>
    class Computation
    {
      const Pde& _eq;
      public:
      Computation(const Pde& eq) : _eq{eq} {}

      Mat<n_state> state;
      void fetch_state(int stride, const double* data)
      {
        for (int i_var = 0; i_var < n_dim + 2; ++i_var) state(i_var) = data[i_var*stride];
        state(n_dim + 2) = data[(2*(n_dim + 2) + 1)*stride];
      }
      Mat<n_update> update_state;
      void fetch_extrap_state(int stride, const double* data)
      {
        for (int i_var = 0; i_var < n_dim + 2; ++i_var) state(i_var) = data[i_var*stride];
        state(n_dim + 2) = 0.;
        update_state = state(Eigen::seqN(0, n_update));
      }

      double mass;
      double kin_ener;
      double pressure;
      void compute_scalars_conv()
      {
        mass = state(n_dim);
        kin_ener = 0;
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) kin_ener += state(i_dim)*state(i_dim);
        kin_ener *= .5/mass;
        pressure = (heat_rat - 1.)*(state((n_dim + 1)) - kin_ener);
      }

      Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
      Mat<n_update, n_dim_flux> flux_conv;
      void compute_flux_conv()
      {
        compute_scalars_conv();
        for (int i_dim = 0; i_dim < n_dim_flux; ++i_dim) {
          flux_conv(n_dim, i_dim) = 0;
          for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
            flux_conv(n_dim, i_dim) += state(j_dim)*normal(j_dim, i_dim);
          }
          double vol_flux = flux_conv(n_dim, i_dim)/mass;
          flux_conv(n_dim + 1, i_dim) = (state((n_dim + 1)) + pressure)*vol_flux;
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
      void compute_scalars_diff()
      {
        bulk_av = std::abs(state(n_dim + 2));
        laplacian_av = 0;
        sqrt_temp = std::sqrt(std::max((state(n_dim + 1) - kin_ener)/mass, 0.)*(heat_rat - 1)/constants::specific_gas_air);
        dyn_visc_coef = _eq.dyn_visc.coefficient(sqrt_temp);
        therm_cond_coef = _eq.therm_cond.coefficient(sqrt_temp);
        energy_cond = therm_cond_coef*(heat_rat - 1)/constants::specific_gas_air;
      }

      Mat<n_extrap, n_dim> gradient;
      Mat<n_update, n_dim> flux_diff;
      void compute_flux_diff()
      {
        if constexpr (n_dim_flux == n_dim) {
          compute_scalars_diff();
          auto seq = Eigen::seqN(0, n_dim);
          auto mmtm = state(seq);
          Mat<n_dim> veloc = mmtm/mass;
          Mat<n_dim, n_dim> veloc_grad = (gradient(seq, all) - veloc*gradient(n_dim, all))/mass;
          Mat<n_dim, n_dim> stress = dyn_visc_coef*(veloc_grad + veloc_grad.transpose())
                                     + ((!_eq._laplacian)*bulk_av*mass - 2./3.*dyn_visc_coef)*veloc_grad.trace()*Mat<n_dim, n_dim>::Identity();
          flux_diff(seq, all) = -stress;
          flux_diff(n_dim, all).setZero();
          Mat<1, n_dim> int_ener_grad = -state(n_dim + 1)/mass/mass*gradient(n_dim, all) + gradient(n_dim + 1, all)/mass - veloc.transpose()*veloc_grad;
          flux_diff(n_dim + 1, all) = -veloc.transpose()*stress - energy_cond*int_ener_grad;
          if (_eq._laplacian) flux_diff -= bulk_av*gradient;
          flux_diff = flux_diff*normal;
        } else HEXED_ASSERT(false, "`compute_flux_diff` requires `n_dim == n_dim_flux`");
      }

      double char_speed;
      void compute_char_speed()
      {
        const double sound_speed = std::sqrt(heat_rat*(heat_rat - 1)*state(n_dim + 1)/state(n_dim)); // numerical estimate (not less than actual speed of sound)
        const double speed = state(Eigen::seqN(0, n_dim)).norm()/state(n_dim);
        char_speed = sound_speed + speed;
      }

      double diffusivity;
      void compute_diffusivity()
      {
        compute_scalars_conv();
        compute_scalars_diff();
        diffusivity = std::abs(laplacian_av) + std::max(std::abs(bulk_av) + dyn_visc_coef/mass, energy_cond/mass);
      }
    };

    /*!
     * Decomposes state vectors into characteristics
     * which are eigenvectors of the Jacobian of the inviscid flux function.
     * This is useful for characteristic-based boundary conditions.
     */
    class Characteristics
    {
      static constexpr int n_var = n_dim + 2;
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
      Characteristics(Mat<n_var> state, Mat<n_dim> direction)
      : dir{direction/direction.norm()},
        mass{state(n_dim)},
        veloc{state(Eigen::seqN(0, n_dim))/mass}
      {
        // compute more properties of the reference state
        double vsq = veloc.squaredNorm();
        double pres = .4*(state(n_dim + 1) - .5*mass*vsq);
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
          state(n_dim + 1) - veloc.dot(mmtm_correction); //! \todo shouldn't this have a `0.5*`?
        // decompose 1D state into eigenvectors
        Mat<1, 3> eig_basis = fact.solve(state_1d).transpose();
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

/*!
 * represents the nonuniform linear advection equation
 * used in the smoothness-based artificial viscosity scheme
 */
template <int n_dim>
class Advection
{
  public:
  static constexpr bool has_diffusion = false;
  static constexpr bool has_convection = true;
  static constexpr bool has_source = true;
  static constexpr int n_state = n_dim + 1;
  static constexpr int n_extrap = n_dim + 1;
  static constexpr int n_update = 1;
  const int _i_node;
  const double _advect_length;
  const double _node;
  const int _adv_var;

  Advection(int i_node, double advect_length, double node)
  : _i_node{i_node}, _advect_length{advect_length}, _node{node}, _adv_var{2*(n_dim + 2) + 3 + 4 + _i_node}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const
  {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_dim; ++i_var) extrap(i_var) = data[i_var*stride];
    extrap(n_dim) = data[_adv_var*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool is_critical) const
  {
    double pseudo = 1 + data[2*(n_dim + 2)*stride]*2/_advect_length;
    double& d = data[_adv_var*stride];
    if (is_critical) d = (d + update(0))/pseudo;
    else d += update(0)/pseudo;
  }

  template <int n_dim_flux>
  class Computation
  {
    const Advection& _eq;
    public:
    Computation(const Advection& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data)
    {
      state = _eq.fetch_extrap(stride, data);
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data)
    {
      for (int i_var = 0; i_var < n_dim + 1; ++i_var) state(i_var) = data[i_var*stride];
      update_state(0) = state(n_dim);
    }

    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_update, n_dim_flux> flux_conv;
    void compute_flux_conv()
    {
      for (int i_dim = 0; i_dim < n_dim_flux; ++i_dim) {
        double nrml_veloc = 0;
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          nrml_veloc += state(j_dim)*normal(j_dim, i_dim);
        }
        flux_conv(0, i_dim) = nrml_veloc*state(n_dim);
      }
    }

    Mat<n_update> source;
    void compute_source()
    {
      source(0) = 2/(_eq._advect_length*_eq._node);
    }

    double char_speed;
    void compute_char_speed()
    {
      char_speed = std::max(1., state(Eigen::seqN(0, n_dim)).norm());
    }

  };
};

/*!
 * represents the uniform linear diffusion equation
 * used in the smoothness-based artificial viscosity scheme
 */
template <int n_dim>
class Smooth_art_visc
{
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = false;
  static constexpr bool has_source = true;
  static constexpr int n_state = 2;
  static constexpr int n_extrap = 1;
  static constexpr int n_update = 1;
  const int _i_step;
  const int _forcing_start;
  const double _diff_time;
  const double _cheby;

  Smooth_art_visc(int i_step, double diff_time, double chebyshev_step)
  : _i_step{i_step}, _forcing_start{2*(n_dim + 2) + 3}, _diff_time{diff_time}, _cheby{chebyshev_step}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const
  {
    Mat<n_extrap> extrap;
    extrap(0) = data[(_forcing_start + _i_step + 1)*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const
  {
    double& d = data[(_forcing_start + _i_step + 1)*stride];
    d += update(0);
    if (critical) d /= 1 + data[2*(n_dim + 2)*stride]*_cheby/_diff_time;
  }

  template <int n_dim_flux>
  class Computation
  {
    const Smooth_art_visc& _eq;
    public:
    Computation(const Smooth_art_visc& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data)
    {
      int f = (_eq._forcing_start + _eq._i_step)*stride;
      state(0) = data[f + stride];
      state(1) = data[f];
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data)
    {
      update_state(0) = *data;
    }

    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_extrap, n_dim> gradient;
    Mat<n_update, n_dim> flux_diff;
    void compute_flux_diff()
    {
      if constexpr (n_dim_flux == n_dim) flux_diff.noalias() = -gradient*normal;
      else HEXED_ASSERT(false, "`compute_flux_diff` requires `n_dim == n_dim_flux`");
    }

    double diffusivity;
    void compute_diffusivity()
    {
      diffusivity = 1;
    }

    Mat<n_update> source;
    void compute_source()
    {
      double f = std::abs(state(1));
      source(0) = ((_eq._i_step == 1) ? std::sqrt(f) : f)/_eq._diff_time;
    }
  };
};

/*!
 * represents the uniform linear diffusion equation
 * used for fixing thermodynamic admissibility
 */
template <int n_dim>
class Fix_therm_admis
{
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = false;
  static constexpr bool has_source = false;
  static constexpr int n_state = n_dim + 2;
  static constexpr int n_update = n_state;
  static constexpr int n_extrap = n_state;

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const
  {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_extrap; ++i_var) extrap(i_var) = data[i_var*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const
  {
    for (int i_var = 0; i_var < n_update; ++i_var) data[i_var*stride] += update(i_var);
  }

  template <int n_dim_flux>
  class Computation
  {
    const Fix_therm_admis& _eq;
    public:
    Computation(const Fix_therm_admis& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data)
    {
      for (int i_var = 0; i_var < n_state; ++i_var) state(i_var) = data[i_var*stride];
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data)
    {
      for (int i_var = 0; i_var < n_state; ++i_var) update_state(i_var) = data[i_var*stride];
    }

    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_extrap, n_dim> gradient;
    Mat<n_update, n_dim> flux_diff;
    void compute_flux_diff()
    {
      if constexpr (n_dim_flux == n_dim) flux_diff.noalias() = -gradient*normal;
      else HEXED_ASSERT(false, "`compute_flux_diff` requires `n_dim == n_dim_flux`");
    }

    double diffusivity;
    void compute_diffusivity()
    {
      diffusivity = 1;
    }
  };
};

}
#endif
