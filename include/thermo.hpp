#ifndef HEXED_THERMO_HPP_
#define HEXED_THERMO_HPP_

#include <cmath>
#include "assert.hpp"

namespace hexed::thermo
{

bool admissible(const double* data, const int n_dim, const int n_qpoint, double& min_mass, double& max_mass, double heat_rat = 1.4);

}

#endif
