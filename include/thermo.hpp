#ifndef HEXED_THERMO_HPP_
#define HEXED_THERMO_HPP_

#include <cmath>
#include "assert.hpp"

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
  HEXED_ASSERT(mass > 0., "nonpositive density") \
  HEXED_ASSERT(ener >= 0., "negative energy") \
  HEXED_ASSERT(pres >= 0., "negative pressure") \

#endif
