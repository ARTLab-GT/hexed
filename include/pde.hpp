#ifndef HEXED_PDE_HPP_
#define HEXED_PDE_HPP_

#include "math.hpp"
#include "Storage_params.hpp"
#include "Gauss_legendre.hpp"

/*! \brief This namespace contains classes representing the different PDEs Hexed can solve.
 * \details They are all possible arguments to the `Spatial` class template.
 * They define the organization of the state data, fluxes, and speeds of information
 * propagation for computing time steps.
 */
namespace hexed::pde {

/*!
 * represents the nonuniform linear advection equation
 * used in the smoothness-based artificial viscosity scheme
 */
template <int n_dim, int row_size>
class Advection {
  const int _n_var;
  static constexpr int _n_adv = row_size;
  Mat<_n_adv> _nodes;
  int _offset;

  public:
  static constexpr bool has_diffusion = false;
  static constexpr bool has_convection = true;
  static constexpr bool has_source = true;
  static constexpr int n_state = n_dim + _n_adv + 1;
  static constexpr int n_extrap = n_dim + _n_adv;
  static constexpr int n_update = _n_adv;
  static constexpr double regular_scale = .1;

  Advection(int n_var, double advect_length, int offset)
  : _n_var{n_var}
  , _offset{offset}
  , _nodes{2*Gauss_legendre(_n_adv).nodes() + Mat<_n_adv>::Constant(math::sign(offset)*.707/_n_adv - 1.)}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_dim; ++i_var) extrap(i_var) = data[i_var*stride];
    for (int i_adv = 0; i_adv < _n_adv; ++i_adv) {
      extrap(n_dim + i_adv) = data[(Storage_params::advection_offset(_n_var) + _offset*_n_adv + i_adv)*stride];
    }
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool is_critical) const {
    double pseudo = 1 + data[Storage_params::tss_offset(_n_var)*stride]*2/data[Storage_params::laplacian_av_offset(_n_var)*stride];
    for (int i_adv = 0; i_adv < _n_adv; ++i_adv) {
      double& d = data[(Storage_params::advection_offset(_n_var) + _offset*_n_adv + i_adv)*stride];
      if (is_critical) d = (d + update(i_adv))/pseudo;
      else d += update(i_adv)/pseudo;
    }
  }

  template <int n_dim_flux>
  class Computation {
    const Advection& _eq;
    public:
    Mat<config::debug_variables> debug_variables; //!< \brief can be populated at any time
    bool debug_vars_set = false;
    Computation(const Advection& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_dim; ++i_var) state(i_var) = data[i_var*stride];
      for (int i_adv = 0; i_adv < _n_adv; ++i_adv) {
        state(n_dim + i_adv) = data[(Storage_params::advection_offset(_eq._n_var) + _eq._offset*_n_adv + i_adv)*stride];
      }
      state(n_dim + _n_adv) = data[Storage_params::laplacian_av_offset(_eq._n_var)*stride];
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_extrap; ++i_var) state(i_var) = data[i_var*stride];
      state(n_extrap) = 0.;
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
      double l = state(n_dim + _n_adv);
      for (int i_adv = 0; i_adv < _n_adv; ++i_adv) {
        double s = state(n_dim + i_adv) - 1.;
        source(i_adv) = 2/l*(1. - math::pow(s*l/regular_scale, 5));
      }
    }

    double char_speed;
    void compute_char_speed() {
      char_speed = (1. + .707/_n_adv)*std::max(1., state(Eigen::seqN(0, n_dim)).norm());
    }

    double decay;
    void compute_decay() {
      double l = state(n_dim + _n_adv);
      decay = 0;
      for (int i_adv = 0; i_adv < _n_adv; ++i_adv) {
        double s = state(n_dim + i_adv) - 1.;
        decay = std::max(decay, 30*2/l*5*l/regular_scale*math::pow(s*l/regular_scale, 4));
      }
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
  static constexpr int n_state = 4 + 1;
  static constexpr int n_extrap = 3;
  static constexpr int n_update = 3;
  const double _diff_time;
  const double _cheby;

  Smooth_art_visc(int n_var, double diff_time, double chebyshev_step)
  : _n_var{n_var}, _diff_time{diff_time}, _cheby{chebyshev_step}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_extrap; ++i_var) extrap(i_var) = data[(Storage_params::forcing_offset(_n_var) + 1 + i_var)*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
    double pseudo = 1 + data[Storage_params::tss_offset(_n_var)*stride]*_cheby
                        /(_diff_time*math::pow(data[Storage_params::laplacian_av_offset(_n_var)*stride], 2));
    for (int i_var = 0; i_var < n_update; ++i_var) {
      double& d = data[(Storage_params::forcing_offset(_n_var) + 1 + i_var)*stride];
      d += update(i_var);
      if (critical) d /= pseudo;
    }
  }

  template <int n_dim_flux>
  class Computation {
    const Smooth_art_visc& _eq;
    public:
    Mat<config::debug_variables> debug_variables; //!< \brief can be populated at any time
    bool debug_vars_set = false;
    Computation(const Smooth_art_visc& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {
      for (int i_var = 0; i_var < n_extrap; ++i_var) state(i_var) = data[(Storage_params::forcing_offset(_eq._n_var) + i_var)*stride];
      state(4) = data[Storage_params::laplacian_av_offset(_eq._n_var)*stride];
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
        source(i_var) = ((i_var == 1) ? std::sqrt(f) : f)/(_eq._diff_time*state(4)*state(4));
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
template <int n_scalar>
class Fix_nonphysical {
  public:
  template <int n_dim, int row_size>
  class Pde {
    public:
    static constexpr bool has_diffusion = true;
    static constexpr bool has_convection = false;
    static constexpr bool has_source = false;
    static constexpr int n_state = n_dim + n_scalar;
    static constexpr int n_update = n_state;
    static constexpr int n_extrap = n_state;

    Pde(int n_var) {}

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
      const Pde& _eq;
      public:
      Mat<config::debug_variables> debug_variables; //!< \brief can be populated at any time
      bool debug_vars_set = false;
      Computation(const Pde& eq) : _eq{eq} {}

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
};

template <int n_dim, int row_size>
class Eikonal {
  int _n_var;
  double _smoothing;
  double _grad_smoothing;
  double _base_diff;
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = true;
  static constexpr bool has_source = true;
  static constexpr int n_state = n_dim + 2;
  static constexpr int n_update = 1;
  static constexpr int n_extrap = n_dim + 1;

  Eikonal(int n_var, double smoothing, double grad_smoothing, double base_diffusion)
  : _n_var{n_var}
  , _smoothing{smoothing}
  , _grad_smoothing{grad_smoothing}
  , _base_diff{base_diffusion}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
    Mat<n_extrap> extrap;
    extrap(0) = data[Storage_params::laplacian_av_offset(_n_var)*stride];
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      extrap(1 + i_dim) = data[(Storage_params::residual_cache_offset(_n_var, row_size) + 1 + i_dim)*stride];
    }
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
    data[Storage_params::laplacian_av_offset(_n_var)*stride] += update(0);
  }

  template <int n_dim_flux>
  class Computation {
    const Eikonal& _eq;
    public:
    Mat<config::debug_variables> debug_variables;
    bool debug_vars_set = false;

    Computation(const Eikonal& eq) : _eq{eq} {}

    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {
      state(0) = data[Storage_params::laplacian_av_offset(_eq._n_var)*stride];
      for (int i_var = 1; i_var < n_state; ++i_var) {
        state(i_var) = data[(Storage_params::residual_cache_offset(_eq._n_var, row_size) + i_var)*stride];
      }
    }
    Mat<n_update> update_state;
    void fetch_extrap_state(int stride, const double* data) {
      update_state(0) = data[0];
      for (int i_var = 0; i_var < n_extrap; ++i_var) state(i_var) = data[i_var*stride];
      state(n_dim + 1) = 0;
    }

    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_update, n_dim_flux> flux_conv;
    void compute_flux_conv() {flux_conv = state(Eigen::seqN(1, n_dim)).transpose()*normal*state(0);}
    Mat<n_extrap, n_dim> gradient;
    Mat<n_update, n_dim_flux> flux_diff;
    void compute_flux_diff() {
      flux_diff.setZero();
      for (int i_dim = 0; i_dim < n_dim_flux; ++i_dim) {
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          flux_diff(0, i_dim) += _eq._smoothing*(state(1 + j_dim) - gradient(0, j_dim))*normal(j_dim, i_dim);
        }
      }
    }
    double char_speed;
    void compute_char_speed() {
      char_speed = 0;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) char_speed += state(1 + i_dim)*state(1 + i_dim);
      char_speed = std::max(1., std::sqrt(char_speed));
    }
    double diffusivity;
    void compute_diffusivity() {
      double grad_sq = 0;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) grad_sq += state(1 + i_dim)*state(1 + i_dim);
      double extra_diff = math::pow(std::sqrt(grad_sq) - 1., 2);
      diffusivity = (1. + _eq._smoothing)*std::abs(state(0)) + _eq._base_diff + _eq._grad_smoothing*extra_diff + _eq._smoothing;
    }

    Mat<n_update> source;
    void compute_source() {
      double grad_sq = 0;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) grad_sq += state(1 + i_dim)*state(1 + i_dim);
      double extra_diff = math::pow(std::sqrt(grad_sq) - 1., 2);
      source(0) = 1. + ((1. + _eq._smoothing)*std::abs(state(0)) + _eq._grad_smoothing*extra_diff)*state(n_dim + 1);
    }
    double decay;
    void compute_decay() {decay = 0.;}
  };
};

template <int n_dim, int row_size>
class Laplace {
  int _read_offset;
  int _write_offset;
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = false;
  static constexpr bool has_source = false;
  static constexpr int n_state = 1;
  static constexpr int n_update = 1;
  static constexpr int n_extrap = 1;

  Laplace(int n_var, int read_offset, int write_offset) : _read_offset{read_offset}, _write_offset{write_offset} {}
  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {return Mat<1>{data[_read_offset*stride]};}
  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
    data[_write_offset*stride] = critical*data[_write_offset*stride] + update(0);
  }

  template <int n_dim_flux>
  class Computation {
    const Laplace& _eq;
    public:
    Mat<config::debug_variables> debug_variables;
    bool debug_vars_set = false;
    Computation(const Laplace& eq) : _eq{eq} {}
    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {state(0) = data[_eq._read_offset*stride];}
    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_extrap, n_dim> gradient;
    Mat<n_update, n_dim_flux> flux_diff;
    void compute_flux_diff() {flux_diff.noalias() = -gradient*normal;}
    double diffusivity;
    void compute_diffusivity() {diffusivity = 1.;}
  };
};

template <int n_dim, int row_size>
class Gradient {
  int _read_offset;
  int _write_offset;
  public:
  static constexpr bool has_diffusion = true; // just to make the spatial kernel compute the gradient
  static constexpr bool has_convection = false;
  static constexpr bool has_source = true;
  static constexpr int n_state = 1;
  static constexpr int n_update = n_dim;
  static constexpr int n_extrap = 1;

  Gradient(int n_var, int read_offset, int write_offset) : _read_offset{read_offset}, _write_offset{write_offset} {}
  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {return Mat<1>{data[_read_offset*stride]};}
  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
    for (int i_var = 0; i_var < n_update; ++i_var) {
      double& d = data[(_write_offset + i_var)*stride];
      d = critical*d + update(i_var);
    }
  }

  template <int n_dim_flux>
  class Computation {
    const Gradient& _eq;
    public:
    Mat<config::debug_variables> debug_variables;
    bool debug_vars_set = false;
    Computation(const Gradient& eq) : _eq{eq} {}
    Mat<n_state> state;
    void fetch_state(int stride, const double* data) {state(0) = data[_eq._read_offset*stride];}
    Mat<n_dim, n_dim_flux> normal = Mat<n_dim, n_dim_flux>::Identity();
    Mat<n_extrap, n_dim> gradient;
    Mat<n_update, n_dim_flux> flux_diff;
    void compute_flux_diff() {flux_diff.setZero();}
    Mat<n_update> source;
    void compute_source() {
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) source(i_dim) = gradient(0, i_dim);
    }
    double diffusivity;
    void compute_diffusivity() {diffusivity = 0.;}
    double decay;
    void compute_decay() {decay = 1.;}
  };
};

}
#endif
