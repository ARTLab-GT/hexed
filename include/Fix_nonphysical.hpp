#ifndef HEXED_FIX_NONPHYSICAL_HPP_
#define HEXED_FIX_NONPHYSICAL_HPP_

#include "math.hpp"
#include "Storage_params.hpp"

namespace hexed {

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

}
#endif
