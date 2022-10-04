#ifndef HEXED_THERMO_HPP_
#define HEXED_THERMO_HPP_

#include <algorithm>
#include <cmath>

namespace hexed::thermo
{

const double mass_ref = 1.225;
const double ener_ref = 101325/.4;
const double min_sound_speed = 1.;

}

#define HEXED_COMPUTE_SCALARS \
  double pres = 0; \
  double mmtm_sq = 0.; \
  for (int j_dim = 0; j_dim < n_var - 2; ++j_dim) { \
    mmtm_sq += READ(j_dim)*READ(j_dim); \
  } \
  pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*mmtm_sq/READ(n_var - 2)); \
  double mass = READ(n_var - 2); \
  double ener = READ(n_var - 1); \
  if (!((pres > 0) && (mass > 0.) && (ener > 0.))) { /* not equiv. to `pres <= 0` cause NaN */ \
    pres = 1e-12; \
    mass = std::max(1e-12*thermo::mass_ref, std::max(mass, std::sqrt(thermo::mass_ref/thermo::ener_ref*.5*mmtm_sq))); \
    ener = std::max(1e-12*thermo::ener_ref, std::max(ener, std::sqrt(thermo::ener_ref/thermo::mass_ref*.5*mmtm_sq))); \
  } \
  { \
    auto loc = " in `" + std::string(__func__) + "`!"; \
    if (mass <= 0.) throw std::runtime_error("nonpositive density" + loc); \
    if (ener <= 0.) throw std::runtime_error("nonpositive energy" + loc); \
    if (!(pres >= 0.)) throw std::runtime_error("nonpositive pressure" + loc); \
  } \

#endif
