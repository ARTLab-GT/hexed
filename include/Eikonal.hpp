#ifndef HEXED_EIKONAL_HPP_
#define HEXED_EIKONAL_HPP_

#include "math.hpp"
#include "Storage_params.hpp"

namespace hexed {

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
    void compute_flux_diff() {flux_diff.setZero();}
    double char_speed;
    void compute_char_speed() {
      char_speed = 0;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) char_speed += state(1 + i_dim)*state(1 + i_dim);
      char_speed = std::max(1., std::sqrt(char_speed));
    }
    double diffusivity;
    void compute_diffusivity() {
      diffusivity = _eq._grad_smoothing + (1. + _eq._smoothing)*std::abs(state(0)) + _eq._base_diff + _eq._smoothing;
    }

    Mat<n_update> source;
    void compute_source() {
      double grad_diff = std::min(1., (state(Eigen::seqN(1, n_dim)) - gradient(0, all).transpose()).norm());
      source(0) = 1. + ((1. + _eq._smoothing)*std::abs(state(0)) + _eq._grad_smoothing*grad_diff)*state(n_dim + 1);
    }
    double decay;
    void compute_decay() {decay = 0.;}
  };
};

}
#endif
