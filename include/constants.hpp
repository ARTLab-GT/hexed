#ifndef HEXED_CONSTANTS_HPP_
#define HEXED_CONSTANTS_HPP_

/*! \file constants.hpp
 * Definition of mathematical and physical constants as `const` global variables
 * (in terms of [SI base units](https://en.wikipedia.org/wiki/SI_base_unit) when possible).
 * Note the doc descriptions with the mathematical symbols for the variables are mostly just
 * so the variables can be linked to in Doxygen but maybe they will be useful to someone...
 */

#include "math.hpp"

namespace hexed
{

//! \name mathematical constants
//!\{
const double pi = M_PI; //!< \brief \f$ \pi \f$
const double degree = pi/360;//!< \brief \f$ ^{\circ} \f$
//!\}

//! \name SI definitions (wikipedia)
//!\{
const double light_speed = 299792458; //!< \brief \f$ c \f$
const double plank = 6.62607015e-34; //!< \brief \f$ h \f$
const double boltzmann = 1.380649e-23; //!< \brief \f$ k \f$
const double avogadro = 6.02214076e23; //!< \brief \f$ mol \f$
//!\}

//! \name empirical values
//!\{
const double mol_mass_air = 28.9647e-3; //!< \brief \f$ \mu_{air} \f$ (from engineering toolbox)
//!\}

//! \name derived quantities
//!\{
#define P(base, exp) math::pow((base), (exp))
const double universal_gas = boltzmann*avogadro; //!< \brief \f$ R_u \f$
const double specific_gas_air = universal_gas/mol_mass_air; //!< \brief \f$ R_{air} \f$
const double stefan_boltzmann = 2*P(pi, 5)*P(boltzmann, 4)/(15*P(light_speed, 2)*P(plank, 3)); //!< \brief \f$ \sigma \f$
#undef P
//!\}

}
#endif
