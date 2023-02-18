#ifndef HEXED_CONSTANTS_HPP_
#define HEXED_CONSTANTS_HPP_

/*! \file constants.hpp
 * Definition of physical constants as `const` global variables
 * (in terms of [SI base units](https://en.wikipedia.org/wiki/SI_base_unit) when possible).
 */

#include "math.hpp"

namespace hexed
{

// SI definitions (wikipedia)
const double light_speed = 299792458;
const double plank = 6.62607015e-34;
const double boltzmann = 1.380649e-23;
const double avogadro = 6.02214076e23;

// empirical values
const double mol_mass_air = 28.9647e-3; // engineering toolbox

// derived quantities
#define P(base, exp) math::pow((base), (exp))
const double universal_gas = boltzmann*avogadro;
const double specific_gas_air = universal_gas/mol_mass_air;
const double stefan_boltzmann = 2*P(M_PI, 5)*P(boltzmann, 4)/(15*P(light_speed, 2)*P(plank, 3));
#undef P

}
#endif
