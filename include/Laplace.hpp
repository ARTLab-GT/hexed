#ifndef HEXED_LAPLACE_HPP_
#define HEXED_LAPLACE_HPP_

#include "math.hpp"
#include "Storage_params.hpp"

namespace hexed {

template <int n_dim, int row_size>
class Laplace {
  int _read_offset;
  int _write_offset;
  public:
  static constexpr bool has_diffusion = true;
  static constexpr bool has_convection = false;
  static constexpr bool has_source = false;
  static constexpr bool needs_size = false;
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

}
#endif
