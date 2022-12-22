#ifndef HEXED_PDE_HPP_
#define HEXED_PDE_HPP_

#include "math.hpp"

namespace hexed::pde
{

#define HEXED_ASSERT_THERM_ADMIS \
  HEXED_ASSERT(state(n_dim) > 0, "nonpositive density"); \
  HEXED_ASSERT(state(n_dim + 1) >= 0, "negative energy"); \
  HEXED_ASSERT(pres >= 0, "negative pressure"); \

template <int n_dim>
class Euler
{
  public:
  Euler() = delete;
  static constexpr int n_var = n_dim + 2;
  static constexpr int curr_start = 0;
  static constexpr int ref_start = n_var;
  static constexpr int n_update = n_var;
  static constexpr double heat_rat = 1.4;

  static constexpr double pressure(Mat<n_var> state)
  {
    double mmtm_sq = 0.;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      mmtm_sq += (state(j_dim))*(state(j_dim));
    }
    return (heat_rat - 1.)*((state(n_dim + 1)) - 0.5*mmtm_sq/(state(n_dim)));
  }

  static constexpr Mat<n_var> flux(Mat<n_var> state, Mat<n_dim> normal)
  {
    Mat<n_var> f;
    f(n_dim) = 0.;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      f(n_dim) += state(j_dim)*normal(j_dim);
    }
    double scaled = f(n_dim)/state(n_dim);
    double pres = pressure(state);
    HEXED_ASSERT_THERM_ADMIS
    f(n_var - 1) = (state(n_dim + 1) + pres)*scaled;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      f(j_dim) = state(j_dim)*scaled + pres*normal(j_dim);
    }
    return f;
  }
};

template <int n_dim>
class Advection
{
  public:
  Advection() = delete;
  static constexpr int n_var = n_dim + 1;
  static constexpr int curr_start = n_dim;
  static constexpr int ref_start = n_dim + 1;
  static constexpr int n_update = 1;

  static constexpr Mat<1> flux(Mat<n_var> state, Mat<n_dim> normal)
  {
    double nrml_veloc = 0;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      nrml_veloc += state(j_dim)*normal(j_dim);
    }
    return Mat<1>::Constant(nrml_veloc*state(n_dim));
  }
};

}
#endif
