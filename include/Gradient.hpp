#ifndef HEXED_GRADIENT_HPP_
#define HEXED_GRADIENT_HPP_

#include "math.hpp"
#include "Storage_params.hpp"

namespace hexed {

template <int n_dim, int row_size>
class Gradient {
  int _read_offset;
  int _write_offset;
  public:
  static constexpr bool has_diffusion = true; // just to make the spatial kernel compute the gradient
  static constexpr bool has_convection = false;
  static constexpr bool has_source = true;
  static constexpr bool needs_size = false;
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
