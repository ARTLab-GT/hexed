#ifndef HEXED_ADVECTION_HPP_
#define HEXED_ADVECTION_HPP_

#include "math.hpp"
#include "Storage_params.hpp"
#include "Gauss_legendre.hpp"

namespace hexed {

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

}
#endif
