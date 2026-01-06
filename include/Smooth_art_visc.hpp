#ifndef HEXED_SMOOTH_ART_VISC_HPP_
#define HEXED_SMOOTH_ART_VISC_HPP_

#include "math.hpp"
#include "Storage_params.hpp"

namespace hexed {

/*!
 * represents the uniform linear diffusion equation
 * used in the smoothness-based artificial viscosity scheme
 */
template <int n_dim, int row_size>
class Smooth_art_visc {
  const int _n_var;
  const double _diff_time;
  const double _cheby;
  const double _wall_shock_width;
  const double _shock_width_growth;
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = false;
  static constexpr bool has_source = true;
  static constexpr bool needs_size = false;
  static constexpr int n_state = 4 + 1;
  static constexpr int n_extrap = 3;
  static constexpr int n_update = 3;

  Smooth_art_visc(int n_var, double diff_time, double chebyshev_step, double wall_shock_width,
                  double shock_width_growth)
  : _n_var{n_var}
  , _diff_time{diff_time}
  , _cheby{chebyshev_step}
  , _wall_shock_width{wall_shock_width}
  , _shock_width_growth{shock_width_growth}
  {}

  Mat<n_extrap> fetch_extrap(int stride, const double* data) const {
    Mat<n_extrap> extrap;
    for (int i_var = 0; i_var < n_extrap; ++i_var) extrap(i_var) = data[(Storage_params::forcing_offset(_n_var) + 1 + i_var)*stride];
    return extrap;
  }

  void write_update(Mat<n_update> update, int stride, double* data, bool critical) const {
    double l = _wall_shock_width + _shock_width_growth*data[Storage_params::laplacian_av_offset(_n_var)*stride];
    double pseudo = 1 + data[Storage_params::tss_offset(_n_var)*stride]*_cheby/(_diff_time*l*l);
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
      double wall_dist = data[Storage_params::laplacian_av_offset(_eq._n_var)*stride];
      state(4) = _eq._wall_shock_width + _eq._shock_width_growth*wall_dist;
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

}
#endif
