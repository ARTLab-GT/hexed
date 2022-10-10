#ifndef HEXED_THERMO_HPP_
#define HEXED_THERMO_HPP_

#include <algorithm>
#include <cmath>

namespace hexed::thermo
{

bool admissible(const double* data, const int n_dim, const int n_qpoint, double heat_rat = 1.4);

}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#define HEXED_COMPUTE_SCALARS \
  double mmtm_sq = 0.; \
  for (int j_dim = 0; j_dim < n_var - 2; ++j_dim) { \
    mmtm_sq += (READ(j_dim))*(READ(j_dim)); \
  } \
  double pres = (heat_rat - 1.)*((READ(n_var - 1)) - 0.5*mmtm_sq/(READ(n_var - 2))); \
  double mass = (READ(n_var - 2)); \
  double ener = (READ(n_var - 1)); \

#pragma GCC diagnostic pop

#define HEXED_ASSERT_ADMISSIBLE \
  { \
    auto loc = " in `" + std::string(__PRETTY_FUNCTION__) + "`!"; \
    /* note: `!(x > 0.)` is not equivalent to `x <= 0.` because `x` could be NaN */ \
    if (!(mass > 0.)) throw std::runtime_error("nonpositive density" + loc); \
    if (!(ener >= 0.)) throw std::runtime_error("negative energy" + loc); \
    if (!(pres >= 0.)) throw std::runtime_error("negative pressure" + loc); \
  } \

#endif
